"""

"""

# ------------------------------------------------------------------------------
# get_version()


def get_version():
    from importlib import metadata

    try:
        # requires existing installation
        version = metadata.version("hipsta")
    except Exception:
        # fall-back if package is not installed, but run directly
        version = "unknown"

    return version
