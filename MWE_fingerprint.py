# %%
# Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

import FlowCal

# import flowkit as fk


# %%
# Define the path to the data
datapath = "data.fcs"

sample = FlowCal.io.FCSData(datapath)

# %%

channels = sample.channels
ch = sorted(channels)
SSC = "SSC - 773/56 - A1"
FSC = "FSC - 456/51 - D2"

ch_leg = [x[-2:] for x in ch]
CHANNEL_ID = {ch_leg[i]: x for i, x in enumerate(ch)}
CHANNEL_ID["SSC"] = SSC
CHANNEL_ID["FSC"] = FSC
CHANNEL_ID["AR"] = "Aspect Ratio_SSC - 773/56 - A1"
CHANNEL_ID["AR_FSC"] = "Aspect Ratio_FSC - 456/51 - D2"
CHANNEL_ID["C3"] = "488 / 594 - 528/46 - C3"

# %%
data_FCM = sample.copy()
data_FCM = FlowCal.transform.to_rfi(data_FCM)
data_FCM_AR = FlowCal.gate.high_low(
    data_FCM,
    channels=[CHANNEL_ID["AR"], CHANNEL_ID["AR_FSC"]],
    low=0.00000001,
    high=0.9999999,
)


# %%
# Basic gating on the fluorescent channel
THRESH_AOBNOB = 3e3
data_FCM_AR_NOB = FlowCal.gate.high_low(
    data_FCM_AR,
    channels=[CHANNEL_ID["C3"]],
    low=1e2,
    high=THRESH_AOBNOB,
)
data_FCM_AR_AOB = FlowCal.gate.high_low(
    data_FCM_AR,
    channels=[CHANNEL_ID["C3"]],
    low=THRESH_AOBNOB,
    high=1e6,
)

data_dict = {
    "AOB": {"data_FCM_AR": data_FCM_AR_AOB},
    "NOB": {"data_FCM_AR": data_FCM_AR_NOB},
}

MIN_SAMPLE_AMOUNT = 3000  # Set a minimum sample amount for subsampling, normally calulated from the data
np.random.seed(123)

# MIN_SAMPLE_AMOUNT = np.min(
#     [
#         np.min(data_pure_AOB["sample_amount"].values),
#         np.min(data_pure_NOB["sample_amount"].values),
#         np.min(data_coculture["sample_amount_AOB"].values),
#         np.min(data_coculture["sample_amount_NOB"].values),
#     ]
# )


# %%
# Analysis by Heyse_2019 et al.
# and Props_2016 et al.
def get_fingerprinting_channels(elem):
    """
    Extract the fingerprinting channels from the data.
    The fingerprinting channels are:
    - FSC (Forward Scatter)
    - SSC (Side Scatter)
    - C3 (SYTO9 fluorescent channel)
    """
    return np.array(
        elem[:, [CHANNEL_ID["FSC"], CHANNEL_ID["SSC"], CHANNEL_ID["C3"]]].tolist()
    )


def subsample_fingerprint_channels(elem):
    """
    Subsample the fingerprinting channels to a minimum sample amount.
    """

    num_rows_2_sample = np.min(
        [
            elem.shape[0],
            MIN_SAMPLE_AMOUNT,
        ]
    )
    return elem[np.random.choice(elem.shape[0], num_rows_2_sample, replace=False)]


def log_minmax_normalize_fingerprint_channels(elem):
    """
    Log transform and min-max normalize the fingerprinting channels.
    The minimum and maximum values are dependent on the flow cytometer used.
    """
    logtrans = np.log10(elem + 100)  # +100 to avoid log(0)

    a = np.log10(100)
    b = np.log10(1e6)  # +1 to avoid log(0)

    return (logtrans - a) / (b - a)  # min-max normalization to [0, 1] range


def calculate_2D_kde(m1, m2):
    """
    Calculate the 2D KDE for two variables m1 and m2.
    # Set grid to 0-1 range with 128 bins
    # and bandwidth of 0.01 according to Props et al. 2016
    """

    # Grid in 0-1 range 128 bins
    X, Y = np.mgrid[0:1:128j, 0:1:128j]
    positions = np.vstack([X.ravel(), Y.ravel()])

    values = np.vstack([m1, m2])

    kernel = stats.gaussian_kde(values, bw_method=0.01)

    Z = np.reshape(kernel(positions).T, X.shape)

    # fig, ax = plt.subplots()
    # ax.imshow(np.rot90(Z), cmap=plt.cm.Blues,
    #         extent=[0, 1, 0, 1])
    # ax.plot(m1, m2, 'k.', alpha=0.1, markersize=1)
    # ax.set_xlim([0, 1])
    # ax.set_ylim([0, 1])
    # plt.show()

    return Z


def minmax_normaliser(Z):
    """
    Normalize the 2D KDE to the [0, 1] interval.
    """
    Z_min = np.min(Z)
    Z_max = np.max(Z)
    if Z_max - Z_min == 0:
        return Z
    return (Z - Z_min) / (Z_max - Z_min)  # min-max normalization to [0, 1] range


def get_kde_fingerprint(elem, plot=False):
    """
    Get the 2D KDE of the fingerprint channels.
    """

    # Grid in 0-1 range 128 bins
    Z1 = calculate_2D_kde(elem[:, 0], elem[:, 1])  # FSC vs SSC
    Z2 = calculate_2D_kde(elem[:, 0], elem[:, 2])  # FSC vs C3
    Z3 = calculate_2D_kde(elem[:, 1], elem[:, 2])  # SSC vs C3
    # Normalize the KDEs to the the [0, 1] interval
    Z1 = minmax_normaliser(Z1)
    Z2 = minmax_normaliser(Z2)
    Z3 = minmax_normaliser(Z3)
    # Reshape the 2D KDEs to 1D arrays and
    # Concatenate the 1D arrays into a single 1D fingerprint
    Z = np.concatenate((Z1.ravel(), Z2.ravel(), Z3.ravel()))
    assert Z.shape == (3 * 128 * 128,), "Fingerprint shape is not correct"

    # Plot the fingerprint
    if plot == True:
        fig, ax = plt.subplots()
        ax.plot(Z)
        plt.xlabel("Population index [-]")
        plt.ylabel("Phenotypic fingerprint [-]")
        plt.ylim([0, 1])
        plt.show()
    # Return the fingerprint
    return Z


# Calculate the diversity index of the phenotypic fingerprint
def hill_number(signal, q=0):
    """
    Calculate the Hill number for a given signal.
    The Hill number is a measure of diversity that takes into account the abundance of each species.
    It is defined as:
    H(q) = (sum_i p_i^q)^(1/(1-q))
    where p_i is the proportion of species i in the signal and q is the order of the Hill number.
    For q=0, it returns the number of species (species richness).
    For q=1, it returns the exponential of the Shannon entropy (effective number of species).
    For q=2, it returns the inverse Simpson index (dominance).
    """
    # Get the unique values and their counts
    unique, counts = np.unique(signal, return_counts=True)
    # Calculate the proportions
    proportions = counts / np.sum(counts)
    if q == 0:
        return len(unique)  # Species richness
    elif q == 1:
        return np.exp(
            -np.sum(proportions * np.log(proportions + 1e-10))
        )  # Effective number of species
    elif q == 2:
        return 1 / np.sum(proportions**2)  # Inverse Simpson index

    else:
        return (np.sum(proportions**q)) ** (1 / (1 - q))  # General Hill number


def calculate_hill_numbers(elem, orders=[0, 1, 2]):
    """
    Calculate the Hill numbers for a given signal.
    The Hill numbers are calculated for the given orders.
    """
    hill_numbers = []
    for q in orders:
        hill_numbers.append(hill_number(elem, q=q))
    return np.array(hill_numbers)


# %%
# Define the processing steps for the data dictionary using the above functions
processing_steps = [
    ("data_FCM_AR", "fingerprinting_channels", get_fingerprinting_channels),
    (
        "fingerprinting_channels",
        "subsampled_fingerprinting_channels",
        subsample_fingerprint_channels,
    ),
    (
        "subsampled_fingerprinting_channels",
        "log_minmax_normalized_fingerprinting_channels",
        log_minmax_normalize_fingerprint_channels,
    ),
    (
        "log_minmax_normalized_fingerprinting_channels",
        "kde_fingerprint",
        lambda elem: get_kde_fingerprint(elem, plot=True),
    ),
    (
        "kde_fingerprint",
        "hill_numbers",
        calculate_hill_numbers,
    ),
]

for FROM, TO, FUNCTION in processing_steps:
    for _, elem in data_dict.items():
        elem[TO] = FUNCTION(elem[FROM])

# %%
HILL_AOB = data_dict["AOB"][TO]
HILL_NOB = data_dict["NOB"][TO]

print("Hill numbers for AOB population:")
print(HILL_AOB)
print("Hill numbers for NOB population:")
print(HILL_NOB)

# %%
