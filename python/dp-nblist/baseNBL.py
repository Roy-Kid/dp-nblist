from typing import Optional
import numpy as np
from numpy.typing import ArrayLike

from .box import Box


class BaseNBL:

    def reset(self):
        """
        reset the neighbor list means clear up all the results but keep configs
        """
        raise NotImplementedError
    
    def build(self, positions:ArrayLike, box:Optional[Box], rc:Optional[float]):
        """
        build the neighbor list by using new positions and configs

        Args:
            positions (ArrayLike): particle positions
            box (Optional[Box]): geometry of the box
            rc (Optional[float]): radius of cutoff

        Raises:
            NotImplementedError: This method must be implemented in the derived class
        """
        raise NotImplementedError
    
    def update(self, positions:ArrayLike):
        """
        update the neighbor list by using new positions

        Args:
            positions (ArrayLike): particle positions

        Raises:
            NotImplementedError: This method must be implemented in the derived class
        """
        raise NotImplementedError
    
    def get_pairs(self)->np.ndarray:
        """
        get neighbor pairs

        Returns:
            np.ndarray: (N_pairs, 2), where first column is the index of center particle, second column is the index of neighbor particle. First column is always smaller than second column.
        """
        pass