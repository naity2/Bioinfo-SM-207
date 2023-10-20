import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from collections import defaultdict
from scipy.stats import gmean
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test


def calcluate_score(response, tpm_scaled, score_name, gene_names, weights=None):
    """Calculate a score based on the expression of a list of genes."""
    gene_names_list = []
    indices = []
    for i, g in enumerate(list(gene_names)):
        if g in tpm_scaled.columns:
            gene_names_list.append(g)
            indices.append(i)
        else:
            print(f"{g} is missing from the expression data")
    gene_exp = tpm_scaled.loc[:, gene_names_list]

    # average expression of a list of genes
    if weights is not None:
        weights = weights[indices]
        assert len(gene_names_list) == len(weights)

    if len(gene_names_list) > 0:
        response[score_name] = np.average(gene_exp, axis=1, weights=weights)
    else:
        response[score_name] = None

def calcluate_geometric_score(response, tpm, score_name, gene_names, weights=None):
    """Calculate a score based on the expression of a list of genes."""
    gene_names_list = []
    indices = []
    for i, g in enumerate(list(gene_names)):
        if g in tpm.columns:
            gene_names_list.append(g)
            indices.append(i)
        else:
            print(f"{g} is missing from the expression data")
    gene_exp = tpm.loc[:, gene_names_list]

    # average expression of a list of genes
    if weights is not None:
        weights = weights[indices]
        assert len(gene_names_list) == len(weights)

    if len(gene_names_list) > 0:
        response[score_name] = gmean(gene_exp, axis=1, weights=weights)
    else:
        response[score_name] = None


def groupby_median(val, median):
    """Group a value into two groups based on the median unless it is missing."""
    if pd.isnull(val):
        return np.nan
    else:
        return "high" if val > median else "low"


# group based on a column
def group_by_col(response, arms, col, group_col=None):
    """Group a column into two groups based on the median for each arm."""
    if group_col is None:
        group_col = col.replace("score", "group")

    group_list = []
    for arm in arms:
        tmp = response.loc[response["Actual Arm Code"] == arm, :]
        median = tmp[col].median()
        # group by median
        groups = tmp[col].map(lambda x: groupby_median(x, median))
        group_list.append(groups)
    groups = pd.concat(group_list)
    groups = groups.loc[response.index]
    assert all(groups.index == response.index)
    response[group_col] = pd.Categorical(
        groups, ordered=True, categories=["low", "high"]
    )


def plot_survival_simple(
    response,
    params,
    arm,
    sig_name,
    ci_show=True,
    do_logrank_test=True,
    timeline="months",
):
    """Plot survival curves for a given arm and a given signature name that doesn't have good or bad."""
    col = f"{sig_name}_group"
    # edge case: groupby col results in only one group for the given arm
    if len(response.loc[response["Actual Arm Code"] == arm, col].unique()) == 1:
        print(f"Only one group for {arm} - {sig_name}")
        return

    if len(params) > 1:
        fig, axes = plt.subplots(
            1, len(params), figsize=(len(params) * 4, 4), dpi=300, sharex=True, sharey=True
        )
    else:
        fig, ax = plt.subplots(1, 1, figsize=(8, 6), dpi=300)

    fig.suptitle(f"{arm} - {sig_name}")
    for i, param in enumerate(params):
        tmp = response.loc[response["Actual Arm Code"] == arm].copy()
        tmp.dropna(subset=[f"{param} time", f"{param} event"], inplace=True)
        kmf = KaplanMeierFitter()
        if len(params) > 1:
            ax = axes[i]

        ts = []
        es = []
        for name, g in tmp.groupby(col):
            t = g[f"{param} time"]
            e = g[f"{param} event"]
            ts.append(t)
            es.append(e)
            kmf.fit(t, e, label=f"{name.title()}")
            kmf.plot_survival_function(ax=ax, ci_show=ci_show)
        if not (i == len(params) - 1):
            ax.get_legend().remove()
        else:
            sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
        ax.set_title(param)
        # logrank test
        if do_logrank_test:
            results = logrank_test(
                ts[0], ts[1], event_observed_A=es[0], event_observed_B=es[1]
            )
            ax.text(
                0.5,
                0.9,
                f"P = {results.p_value:.2e}",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
        ax.set(xlabel=f"Timeline {timeline}", ylabel="Est. probability")
    plt.tight_layout()


def plot_abundance(data, group, target_col, normalize=True, figsize=(8, 6)):
    """Plot the abundance of each value in the target_col for each group."""
    col_name = "Percent" if normalize else "Count"
    abundance = (
        data.groupby(group)[target_col]
        .value_counts(normalize=normalize)
        .reset_index(name=col_name)
    )
    abundance = abundance.pivot(index=group, columns=target_col, values=col_name)
    # stacked bar plot
    abundance.plot(kind="bar", stacked=True, figsize=figsize).legend(
        loc="center left", bbox_to_anchor=(1, 0.5)
    )
    plt.xticks(rotation=30, ha="right")
    plt.title(target_col.split("_")[0].title())
    sns.despine()


def calculate_indices(response):
    """Calculate the  the B cell and immunosuppression indices from IO360."""
    indices = {
        "BcellProliferationIndex": [
            "Proliferation_IO360",
            "MSI Predictor_IO360",
            "APM Loss_IO360",
            "B Cells_IO360",
            "Glycolytic Activity_IO360",
            "CD45_IO360",
            "JAKSTAT Loss_IO360",
        ],
        "ImmunosuppressionIndex": [
            "MMR Loss_IO360",
            "Hypoxia_IO360",
            "Apoptosis_IO360",
            "NOS2_IO360",
            "MAGEs_IO360",
            "Mast Cells_IO360",
            "TGF-Beta_IO360",
            "ARG1_IO360",
            "Endothelial Cells_IO360",
            "Stroma_IO360",
            "B7-H3_IO360",
            "Myeloid Inflammation_IO360",
        ],
    }

    for i, sigs in indices.items():
        scores = [s + "_score" for s in sigs]
        tmp = response.loc[:, scores].to_numpy()
        rms = np.sqrt(np.mean(tmp**2, axis=1))
        response[i] = rms


def score_group_signatures(response, data_scaled, signatures, arms, score_prefix=""):
    # how many genes are in response
    gene_counts = {}

    sig_names = []

    # P values
    sig_pvals = defaultdict(list)

    # significant signatures
    significant_sigs = defaultdict(list)

    for sig_name, sig_genes in signatures.items():
        gene_names = [
            gene for gene in sig_genes if gene in data_scaled.columns.to_list()
        ]
        gene_counts[sig_name] = len(gene_names)

        # calculate the score if gene_names is not empty
        if gene_names:
            sig_names.append(sig_name)
            calcluate_score(
                response, data_scaled, f"{score_prefix}{sig_name}_score", gene_names
            )
            # group the score into two groups
            group_by_col(response, arms, f"{score_prefix}{sig_name}_score")

            # calculate the logrank test
            for arm in arms:
                tmp = response.loc[response["Actual Arm Code"] == arm].copy()
                # only consider PFS
                tmp.dropna(subset=["PFS time", "PFS event"], inplace=True)
                col = f"{score_prefix}{sig_name}_group"
                ts = []
                es = []
                for _, g in tmp.groupby(col):
                    t = g["PFS time"]
                    e = g["PFS event"]
                    ts.append(t)
                    es.append(e)

                # calculate the logrank test
                res = logrank_test(
                    ts[0], ts[1], event_observed_A=es[0], event_observed_B=es[1]
                )
                sig_pvals[arm].append(res.p_value)
                if res.p_value < 0.05:
                    significant_sigs[arm].append(sig_name)
        else:
            print(f"{sig_name} missing")

    return gene_counts, sig_names, sig_pvals, significant_sigs


def univariate_cox_hr(response, arm, params, formula, num_digits=2):
    cox_res = defaultdict(list)
    for param in params:
        tmp = response.loc[response["Actual Arm Code"] == arm].copy()
        tmp.dropna(subset=[f"{param} time", f"{param} event"], inplace=True)
        cph = CoxPHFitter()
        cph.fit(
            tmp,
            duration_col=f"{param} time",
            event_col=f"{param} event",
            formula=formula,
        )
        hr = cph.summary["exp(coef)"][0]
        lower = cph.summary["exp(coef) lower 95%"][0]
        higher = cph.summary["exp(coef) upper 95%"][0]
        cox_res["Arm"].append(arm)
        cox_res["param"].append(param)
        cox_res["HR"].append(hr)
        cox_res["lower"].append(lower)
        cox_res["higher"].append(higher)

    df = pd.DataFrame(cox_res)
    if num_digits is not None:
        df = df.round(num_digits)

    return df


def logrank_p(response, arm, params, group_col, num_digits=None):
    logrank_results = defaultdict(list)
    for param in params:
        tmp = response.loc[response["Actual Arm Code"] == arm].copy()
        tmp.dropna(subset=[f"{param} time", f"{param} event"], inplace=True)
        ts = []
        es = []
        for _, g in tmp.groupby(group_col):
            t = g[f"{param} time"]
            e = g[f"{param} event"]
            ts.append(t)
            es.append(e)

        # calculate the logrank test
        res = logrank_test(ts[0], ts[1], event_observed_A=es[0], event_observed_B=es[1])
        logrank_results["Arm"].append(arm)
        logrank_results["group_col"].append(group_col)
        logrank_results["param"].append(param)
        logrank_results["P value"].append(res.p_value)

    df = pd.DataFrame(logrank_results)
    if num_digits is not None:
        df = df.round(num_digits)

    return df


def generate_stats(response, arm, params, group_col, num_digits=None):
    # if group contains space, enclose it with ``
    if " " in group_col:
        formula = f"`{group_col}`"
    else:
        formula = group_col
    hr_df = univariate_cox_hr(response, arm, params, formula, num_digits)
    p_df = logrank_p(response, arm, params, group_col, num_digits)

    # merge two dfs
    return pd.merge(p_df, hr_df, on=["Arm", "param"])

def plot_survival_two_arms(response, signame, param="PFS", timeline="years"):
    fig, axes = plt.subplots(1, 2, figsize=(12, 6), dpi=300, sharex=True, sharey=True)
    plt.rc('font', size=20)
    arms = response["Actual Arm Code"].unique().tolist()
    for i, arm in enumerate(arms):
        tmp = response.loc[response["Actual Arm Code"] == arm].copy()
        tmp.dropna(subset = [f"{param} time", f"{param} event"], inplace=True)
        col = f"{signame}_group"
        kmf = KaplanMeierFitter()
        ax = axes[i]
        
        ts = []
        es = []
        for name, g in tmp.groupby(col):
            t = g[f"{param} time"]
            e = g[f"{param} event"]
            ts.append(t)
            es.append(e)
            kmf.fit(t, e, label=f"{name.title()}")
            kmf.plot_survival_function(ax=ax)
        if i==0:
            ax.get_legend().remove()
        else:
            sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
        # logrank test
        results = logrank_test(ts[0], ts[1],
                               event_observed_A=es[0], event_observed_B=es[1])
        ax.text(0.5, 0.9, f"P = {results.p_value:.2e}", 
                ha="center", va="center", transform=ax.transAxes)
        ax.set(xlabel=f'Timeline ({timeline})', ylabel='Est. probability')
        ax.set_title(arm)
    plt.tight_layout()


def plot_survival_ecotyper(
    response,
    params,
    arm,
    sig_name,
    ci_show=False,
    loc = None,
    show_p = False,
    timeline="months",
    xlim = None,
    ylim = None,
):
    col = f"{sig_name}_group"
    # edge case: groupby col results in only one group for the given arm
    if len(response.loc[response["Actual Arm Code"] == arm, col].unique()) == 1:
        print(f"Only one group for {arm} - {sig_name}")
        return

    if len(params) > 1:
        fig, axes = plt.subplots(
            1, len(params), figsize=(len(params) * 4, 4), dpi=300, sharex=True, sharey=True
        )
    else:
        fig, ax = plt.subplots(1, 1, figsize=(8, 6), dpi=300)

    fig.suptitle(f"{arm} - {sig_name}")
    for i, param in enumerate(params):
        tmp = response.loc[response["Actual Arm Code"] == arm].copy()
        tmp.dropna(subset=[f"{param} time", f"{param} event"], inplace=True)
        kmf = KaplanMeierFitter()
        if len(params) > 1:
            ax = axes[i]

        ts = []
        es = []
        num_samples = {}
        for name, g in tmp.groupby(col):
            t = g[f"{param} time"]
            e = g[f"{param} event"]
            num_samples[name] = len(g)
            ts.append(t)
            es.append(e)
            kmf.fit(t, e, label=name)
            kmf.plot_survival_function(ax=ax, ci_show=ci_show, loc=loc)
        ax.set_title(param)        
        ax.set(xlabel=f"Timeline ({timeline})", ylabel="Est. probability")

        # add number of samples for each group in the figure legend
        # Get the legend handles and labels
        handles, labels = ax.get_legend_handles_labels()

        # Modify the labels in a for loop
        for j in range(len(labels)):
            labels[j] = f"{labels[j]} (n={num_samples[labels[j]]})"

        # Add legend to the plot
        if i == len(params) - 1:  
            ax.legend(handles, labels,
                      loc='upper left', bbox_to_anchor=(1, 1))
        else:
            ax.get_legend().remove()

        # logrank test
        if show_p:
            results = logrank_test(ts[0], ts[1],
                                   event_observed_A=es[0],
                                   event_observed_B=es[1])
            ax.text(0.5, 0.9, f"P = {results.p_value:.2e}", 
                    ha="center", va="center", transform=ax.transAxes)

    if xlim is not None:
        plt.xlim(*xlim)
    if ylim is not None:
        plt.ylim(*ylim)
    plt.tight_layout()                                                      