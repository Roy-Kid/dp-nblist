import numpy as np
import math

def hash_function_int(point):
    x, y, z = point
    # 将x、y、z分别转换为32位整数
    x = (x | (x << 16)) & 0x030000FF
    x = (x | (x << 8)) & 0x0300F00F
    x = (x | (x << 4)) & 0x030C30C3
    x = (x | (x << 2)) & 0x09249249

    y = (y | (y << 16)) & 0x030000FF
    y = (y | (y << 8)) & 0x0300F00F
    y = (y | (y << 4)) & 0x030C30C3
    y = (y | (y << 2)) & 0x09249249

    z = (z | (z << 16)) & 0x030000FF
    z = (z | (z << 8)) & 0x0300F00F
    z = (z | (z << 4)) & 0x030C30C3
    z = (z | (z << 2)) & 0x09249249

    # 将x、y、z按位交错合并为一个整数
    return x | (y << 1) | (z << 2)


def hash_function_float_1(point, coefficient):
    x, y, z = point
    # 将浮点数坐标乘以放大倍数，然后取整，得到离散化的整数坐标
    discrete_x = int(x * coefficient)
    discrete_y = int(y * coefficient)
    discrete_z = int(z * coefficient)

    # 对离散化的整数坐标进行位运算，生成哈希值
    x = (discrete_x | (discrete_x << 16)) & 0x030000FF
    x = (x | (x << 8)) & 0x0300F00F
    x = (x | (x << 4)) & 0x030C30C3
    x = (x | (x << 2)) & 0x09249249

    y = (discrete_y | (discrete_y << 16)) & 0x030000FF
    y = (y | (y << 8)) & 0x0300F00F
    y = (y | (y << 4)) & 0x030C30C3
    y = (y | (y << 2)) & 0x09249249

    z = (discrete_z | (discrete_z << 16)) & 0x030000FF
    z = (z | (z << 8)) & 0x0300F00F
    z = (z | (z << 4)) & 0x030C30C3
    z = (z | (z << 2)) & 0x09249249

    # 将离散化后的整数坐标按位交错合并为一个整数，作为最终的哈希值
    hash_value = x | (y << 1) | (z << 2)
    return hash_value

    # # 将离散化后的整数坐标按位交错合并为一个整数，作为最终的哈希值
    # hash_value = x | (y << 1) | (z << 2)

    # hash_bits = int(math.log10(hash_value) + 1) if hash_value != 0 else 1

    # return hash_bits


def hash_function_float_2(point):
    x, y, z = point
    # 将浮点坐标映射到整数坐标范围（例如，0到2^16 - 1）
    max_coord = 1  # 2^16 - 1
    x_int = int(x * max_coord)
    y_int = int(y * max_coord)
    z_int = int(z * max_coord)

    # 将整数坐标合并成一个一维哈希值，使用 Z-order curve 映射
    combined = 0
    for i in range(16):  # 假设每个坐标分量都是 16 位
        combined |= (x_int & 1) << (3 * i)
        combined |= (y_int & 1) << (3 * i + 1)
        combined |= (z_int & 1) << (3 * i + 2)
        x_int >>= 1
        y_int >>= 1
        z_int >>= 1

    return combined


def custom_hash(point, max_bits=32):
    x, y, z = point
    # Ensure x, y, and z are non-negative integers.
    x = int(x)
    y = int(y)
    z = int(z)

    # Calculate the maximum coordinate value based on the number of bits.
    max_coord = 2 ** max_bits - 1

    # Ensure the coordinates are within bounds.
    if x < 0 or x > max_coord or y < 0 or y > max_coord or z < 0 or z > max_coord:
        raise ValueError("Coordinates are out of bounds for the specified number of bits.")

    # Interleave the bits of x, y, and z to create a single integer.
    result = 0
    for i in range(max_bits):
        result |= ((x >> i) & 1) << (3 * i + 2)
        result |= ((y >> i) & 1) << (3 * i + 1)
        result |= ((z >> i) & 1) << (3 * i)

    return result




