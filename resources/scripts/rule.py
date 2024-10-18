"""
Functions for use in Snakemake rule definitions.
"""

import re


def parse_info(info: dict) -> dict:
    """
    Parse dictionary of sample info.

    Arguments:
        ``info``: dictionary of sample info.

    Returns:
        Dictionary of parsed sample info.
    """
    samples = list(info.keys())

    return {"samples": samples}


def get_optional_flags(**kwargs) -> str:
    """
    Get optional flags for a command line script.

    Arguments:
        ``kwargs``: keyword arguments for flags.

    Returns:
        String containing flag to be inserted into shell command.
    """
    assert all(
        not isinstance(value, list) for value in kwargs.values()
    ), "Arguments cannot be lists"
    flags = [
        (
            f"--{str(key).replace('_', '-')}"
            if isinstance(value, bool)
            else (
                f"--{str(key).replace('_', '-')} '{str(value)}'"
                if re.search(r"\s|\"", str(value))
                else f"--{str(key).replace('_', '-')} {str(value)}"
            )
        )
        for key, value in kwargs.items()
        if value
    ]
    return " ".join(flags)


def parse_integration_method(method: str | None) -> str:
    """
    Parse integration method.

    Arguments:
        ``method``: integration method.

    Returns:
        Name of integration function corresponding to chosen method.
    """
    if method is None:
        return None
    match str(method).lower():
        case "rpca":
            return "RPCAIntegration"
        case "jpca":
            return "JPCAIntegration"
        case "cca":
            return "CCAIntegration"
        case "fastmnn":
            return "ReducedMNNIntegration"
        case "harmony":
            return "HarmonyIntegration"
        case "scvi":
            return "scVIIntegration"
        case _:
            raise ValueError(
                f"Integration method {method} not recognized; must be one of 'rpca', 'jpca', 'cca', 'fastmnn', 'harmony', or 'scvi'."
            )
