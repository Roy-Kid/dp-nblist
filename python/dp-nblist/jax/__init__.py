from warnings import WarningMessage


try:
    import torch
except ModuleNotFoundError:
    WarningMessage("jax not installed, skip")