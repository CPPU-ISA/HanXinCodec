# -*- coding: utf-8 -*-
"""汉信码 (Han Xin Code) 编码 / 解码 —— 纯 Python 标准库实现。

标准依据: GB/T 21049-2022《汉信码》/ ISO/IEC 20830:2021 Han Xin Code。
算法: 四角定位图形 + 功能信息(GF(16) RS) + 掩模 + 13 路栅栏交织 + GF(256) RS 纠错。

用法:
    python hanxin.py encode "要编码的文本" -o hanxin.png [--ec L2] [--version 6] [--mask 2]
    python hanxin.py decode hanxin.png [-j]        # -j 输出 JSON
    python hanxin.py decode code.png --ascii       # 同时打印模块矩阵

依赖: 生成只需要标准库; 解析图片需要 Pillow (pip install pillow)。
"""
from __future__ import annotations

import argparse
import json
import os
import struct
import sys
import zlib
from functools import lru_cache

MIN_VERSION = 1
MAX_VERSION = 84
EC_LEVELS = ("L1", "L2", "L3", "L4")


class HanXinError(Exception):
    """汉信码编解码错误基类。"""


class ReedSolomonError(HanXinError):
    """RS 纠错失败(错误个数超出纠错能力)。"""


class InvalidBitStreamError(HanXinError):
    """信息位流不合法。"""


class FunctionInformationError(HanXinError):
    """功能信息不可解码。"""


def dimension_for_version(version: int) -> int:
    if not isinstance(version, int) or not MIN_VERSION <= version <= MAX_VERSION:
        raise ValueError(f"汉信码版本必须是 {MIN_VERSION}..{MAX_VERSION} 之间的整数")
    return 2 * version + 21


def version_for_dimension(dimension: int) -> int:
    if dimension < 23 or dimension > 189 or dimension % 2 == 0:
        raise ValueError("汉信码边长必须是 23..189 之间的奇数")
    return (dimension - 21) // 2


# --------------------------------------------------------------------------
# GF(2^n): 4 次域用于功能信息, 8 次域用于数据区 RS 纠错
# --------------------------------------------------------------------------
class BinaryExtensionField:
    def __init__(self, degree: int, polynomial: int):
        size = 1 << degree
        if not size <= polynomial < 2 * size:
            raise ValueError(f"本原多项式必须是 {degree} 次二进制多项式")
        self.degree = degree
        self.size = size
        self.order = size - 1
        self.polynomial = polynomial
        self.exp_table = [0] * (2 * self.order)
        self.log_table = [-1] * size

        value = 1
        for exponent in range(self.order):
            if self.log_table[value] != -1:
                raise ValueError("多项式对生成元 alpha=2 不是本原的")
            self.exp_table[exponent] = value
            self.log_table[value] = exponent
            value <<= 1
            if value & size:
                value ^= polynomial
            value &= self.order
        if value != 1:
            raise ValueError("多项式不能生成该有限域")
        for exponent in range(self.order, 2 * self.order):
            self.exp_table[exponent] = self.exp_table[exponent - self.order]

    def exp(self, exponent: int) -> int:
        return self.exp_table[exponent % self.order]

    def log(self, value: int) -> int:
        if value == 0:
            raise ValueError(f"GF({self.size}) 中 log(0) 无定义")
        return self.log_table[value]

    def inverse(self, value: int) -> int:
        if value == 0:
            raise ValueError("零元没有乘法逆元")
        return self.exp_table[self.order - self.log_table[value]]

    def multiply(self, left: int, right: int) -> int:
        if left == 0 or right == 0:
            return 0
        return self.exp_table[self.log_table[left] + self.log_table[right]]


GF16 = BinaryExtensionField(4, 0b10011)          # alpha^4 + alpha + 1
GF256 = BinaryExtensionField(8, 0b101100011)     # 数据区 RS 纠错


def _poly_multiply(field: BinaryExtensionField, left, right):
    product = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        if a == 0:
            continue
        for j, b in enumerate(right):
            product[i + j] ^= field.multiply(a, b)
    return product


def rs_generator(field: BinaryExtensionField, degree: int, first_root: int = 1):
    generator = [1]
    for index in range(degree):
        generator = _poly_multiply(field, generator, [1, field.exp(first_root + index)])
    return generator


def rs_encode(field: BinaryExtensionField, data, correction: int, first_root: int = 1):
    data = list(data)
    generator = rs_generator(field, correction, first_root)
    work = data + [0] * correction
    for index in range(len(data)):
        coefficient = work[index]
        if coefficient == 0:
            continue
        for offset in range(1, len(generator)):
            work[index + offset] ^= field.multiply(generator[offset], coefficient)
    return data + work[len(data):]


class _Poly:
    """GF 上的一元多项式, 系数按"最高次在前"存放。"""

    __slots__ = ("field", "coefficients")

    def __init__(self, field: BinaryExtensionField, coefficients):
        coefficients = list(coefficients)
        first = 0
        while first < len(coefficients) - 1 and coefficients[first] == 0:
            first += 1
        self.field = field
        self.coefficients = coefficients[first:]

    @property
    def degree(self) -> int:
        return len(self.coefficients) - 1

    @property
    def is_zero(self) -> bool:
        return self.coefficients[0] == 0

    def coefficient(self, degree: int) -> int:
        return self.coefficients[len(self.coefficients) - 1 - degree]

    def evaluate(self, value: int) -> int:
        if value == 0:
            return self.coefficient(0)
        result = self.coefficients[0]
        for coefficient in self.coefficients[1:]:
            result = self.field.multiply(result, value) ^ coefficient
        return result

    def add(self, other: "_Poly") -> "_Poly":
        if self.is_zero:
            return other
        if other.is_zero:
            return self
        smaller, larger = self.coefficients, other.coefficients
        if len(smaller) > len(larger):
            smaller, larger = larger, smaller
        total = list(larger)
        offset = len(larger) - len(smaller)
        for index, value in enumerate(smaller):
            total[offset + index] ^= value
        return _Poly(self.field, total)

    def multiply(self, other: "_Poly") -> "_Poly":
        if self.is_zero or other.is_zero:
            return _Poly(self.field, [0])
        return _Poly(self.field, _poly_multiply(self.field, self.coefficients, other.coefficients))

    def scale(self, scalar: int) -> "_Poly":
        if scalar == 0:
            return _Poly(self.field, [0])
        if scalar == 1:
            return self
        return _Poly(self.field, [self.field.multiply(c, scalar) for c in self.coefficients])

    def multiply_monomial(self, degree: int, coefficient: int) -> "_Poly":
        if coefficient == 0:
            return _Poly(self.field, [0])
        result = [0] * (len(self.coefficients) + degree)
        for index, value in enumerate(self.coefficients):
            result[index] = self.field.multiply(value, coefficient)
        return _Poly(self.field, result)


def _monomial(field, degree, coefficient):
    coefficients = [0] * (degree + 1)
    coefficients[0] = coefficient
    return _Poly(field, coefficients)


def _euclidean_algorithm(field, left, right, correction):
    a, b = left, right
    if a.degree < b.degree:
        a, b = b, a
    previous_remainder, remainder = a, b
    previous_auxiliary, auxiliary = _Poly(field, [0]), _Poly(field, [1])

    while 2 * remainder.degree >= correction:
        older_remainder, older_auxiliary = previous_remainder, previous_auxiliary
        previous_remainder, previous_auxiliary = remainder, auxiliary
        if previous_remainder.is_zero:
            raise ReedSolomonError("欧几里得算法出现零余式")
        remainder = older_remainder
        quotient = _Poly(field, [0])
        denominator_inverse = field.inverse(previous_remainder.coefficient(previous_remainder.degree))
        while remainder.degree >= previous_remainder.degree and not remainder.is_zero:
            degree_difference = remainder.degree - previous_remainder.degree
            scale = field.multiply(remainder.coefficient(remainder.degree), denominator_inverse)
            quotient = quotient.add(_monomial(field, degree_difference, scale))
            remainder = remainder.add(previous_remainder.multiply_monomial(degree_difference, scale))
        auxiliary = quotient.multiply(previous_auxiliary).add(older_auxiliary)
        if remainder.degree >= previous_remainder.degree:
            raise ReedSolomonError("欧几里得多项式除法未能降次")

    normalization = auxiliary.coefficient(0)
    if normalization == 0:
        raise ReedSolomonError("错误定位多项式无法归一化")
    inverse = field.inverse(normalization)
    return auxiliary.scale(inverse), remainder.scale(inverse)


def _error_locations(field, locator):
    count = locator.degree
    if count == 1:
        return [locator.coefficient(1)]
    locations = []
    for value in range(1, field.size):
        if len(locations) >= count:
            break
        if locator.evaluate(value) == 0:
            locations.append(field.inverse(value))
    if len(locations) != count:
        raise ReedSolomonError("错误定位多项式的根数与次数不符")
    return locations


def _error_magnitudes(field, evaluator, locations, first_root):
    magnitudes = []
    for index, location in enumerate(locations):
        inverse_location = field.inverse(location)
        denominator = 1
        for other, other_location in enumerate(locations):
            if other == index:
                continue
            denominator = field.multiply(
                denominator, field.multiply(other_location, inverse_location) ^ 1
            )
        magnitude = field.multiply(
            evaluator.evaluate(inverse_location), field.inverse(denominator)
        )
        if first_root != 0:
            magnitude = field.multiply(magnitude, inverse_location)
        magnitudes.append(magnitude)
    return magnitudes


def rs_decode(field: BinaryExtensionField, received, correction: int, first_root: int = 1):
    """返回 (纠正后的码字, 纠正的错误码字个数)。"""
    received = list(received)
    if not received or len(received) > field.order:
        raise ReedSolomonError(f"GF({field.size}) 上的码块长度必须在 1..{field.order}")

    polynomial = _Poly(field, received)
    syndromes = [0] * correction
    has_error = False
    for index in range(correction):
        syndrome = polynomial.evaluate(field.exp(index + first_root))
        syndromes[correction - 1 - index] = syndrome
        if syndrome:
            has_error = True
    if not has_error:
        return received, 0

    locator, evaluator = _euclidean_algorithm(
        field, _monomial(field, correction, 1), _Poly(field, syndromes), correction
    )
    locations = _error_locations(field, locator)
    magnitudes = _error_magnitudes(field, evaluator, locations, first_root)
    for location, magnitude in zip(locations, magnitudes):
        position = len(received) - 1 - field.log(location)
        if not 0 <= position < len(received):
            raise ReedSolomonError("算出的错误位置越过码块边界")
        received[position] ^= magnitude

    verified = _Poly(field, received)
    for index in range(correction):
        if verified.evaluate(field.exp(index + first_root)) != 0:
            raise ReedSolomonError("纠错后仍存在非零伴随式")
    return received, len(locations)


# --------------------------------------------------------------------------
# 标准表格 (由 GB/T 21049-2022 表 A.1 / 表 B.1 / 表 18 数值)
# --------------------------------------------------------------------------
# GB/T 21049-2022 表 B.1 码块结构: 版本 -> L1..L4 -> [(块数, 数据码字, 纠错码字)]
# 数值为标准的规范性参数(与 Zint 2.16 backend/hanxin.h 的 Table D1 逐项一致), 独立实现为互操作必须使用同样的数值
RS_BLOCKS = [
    [[(1, 21, 4)], [(1, 17, 8)], [(1, 13, 12)], [(1, 9, 16)]],
    [[(1, 31, 6)], [(1, 25, 12)], [(1, 19, 18)], [(1, 15, 22)]],
    [[(1, 42, 8)], [(1, 34, 16)], [(1, 26, 24)], [(1, 20, 30)]],
    [[(1, 46, 8)], [(1, 38, 16)], [(1, 30, 24)], [(1, 22, 32)]],
    [[(1, 57, 12)], [(1, 49, 20)], [(1, 37, 32)], [(1, 14, 20), (1, 13, 22)]],
    [[(1, 70, 14)], [(1, 58, 26)], [(1, 24, 20), (1, 22, 18)], [(1, 16, 24), (1, 18, 26)]],
    [[(1, 84, 16)], [(1, 70, 30)], [(1, 26, 22), (1, 28, 24)], [(2, 14, 20), (1, 12, 20)]],
    [[(1, 99, 18)], [(1, 40, 18), (1, 41, 18)], [(1, 31, 26), (1, 32, 28)], [(2, 16, 24), (1, 15, 22)]],
    [[(1, 114, 22)], [(2, 48, 20)], [(2, 24, 20), (1, 26, 22)], [(2, 18, 28), (1, 18, 26)]],
    [[(1, 131, 24)], [(1, 52, 22), (1, 57, 24)], [(2, 27, 24), (1, 29, 24)], [(2, 21, 32), (1, 19, 30)]],
    [[(1, 135, 26)], [(1, 56, 24), (1, 57, 24)], [(2, 28, 24), (1, 31, 26)], [(2, 22, 32), (1, 21, 32)]],
    [[(1, 153, 28)], [(1, 62, 26), (1, 65, 28)], [(2, 32, 28), (1, 33, 28)], [(3, 17, 26), (1, 22, 30)]],
    [[(1, 86, 16), (1, 85, 16)], [(1, 71, 30), (1, 72, 30)], [(2, 37, 32), (1, 35, 30)], [(3, 20, 30), (1, 21, 32)]],
    [[(1, 94, 18), (1, 95, 18)], [(2, 51, 22), (1, 55, 24)], [(3, 30, 26), (1, 31, 26)], [(4, 18, 28), (1, 17, 24)]],
    [[(1, 104, 20), (1, 105, 20)], [(2, 57, 24), (1, 61, 26)], [(3, 33, 28), (1, 36, 30)], [(4, 20, 30), (1, 19, 30)]],
    [[(1, 115, 22), (1, 114, 22)], [(2, 65, 28), (1, 61, 26)], [(3, 38, 32), (1, 33, 30)], [(5, 19, 28), (1, 14, 24)]],
    [[(1, 126, 24), (1, 125, 24)], [(2, 70, 30), (1, 69, 30)], [(4, 33, 28), (1, 29, 26)], [(5, 20, 30), (1, 19, 30)]],
    [[(1, 136, 26), (1, 137, 26)], [(3, 56, 24), (1, 59, 26)], [(5, 35, 30)], [(6, 18, 28), (1, 21, 28)]],
    [[(1, 148, 28), (1, 149, 28)], [(3, 61, 26), (1, 64, 28)], [(7, 24, 20), (1, 23, 22)], [(6, 20, 30), (1, 21, 32)]],
    [[(3, 107, 20)], [(3, 65, 28), (1, 72, 30)], [(7, 26, 22), (1, 23, 22)], [(7, 19, 28), (1, 20, 32)]],
    [[(3, 115, 22)], [(4, 56, 24), (1, 63, 28)], [(7, 28, 24), (1, 25, 22)], [(8, 18, 28), (1, 21, 22)]],
    [[(2, 116, 22), (1, 122, 24)], [(4, 56, 24), (1, 72, 30)], [(7, 28, 24), (1, 32, 26)], [(8, 18, 28), (1, 24, 30)]],
    [[(3, 127, 24)], [(5, 51, 22), (1, 62, 26)], [(7, 30, 26), (1, 35, 26)], [(8, 20, 30), (1, 21, 32)]],
    [[(2, 135, 26), (1, 137, 26)], [(5, 56, 24), (1, 59, 26)], [(7, 33, 28), (1, 30, 28)], [(11, 16, 24), (1, 19, 26)]],
    [[(3, 105, 20), (1, 121, 22)], [(5, 61, 26), (1, 57, 26)], [(9, 28, 24), (1, 28, 22)], [(10, 19, 28), (1, 18, 30)]],
    [[(2, 157, 30), (1, 150, 28)], [(5, 65, 28), (1, 61, 26)], [(8, 33, 28), (1, 34, 30)], [(10, 19, 28), (2, 15, 26)]],
    [[(3, 126, 24), (1, 115, 22)], [(7, 51, 22), (1, 54, 22)], [(8, 35, 30), (1, 37, 30)], [(15, 15, 22), (1, 10, 22)]],
    [[(4, 105, 20), (1, 103, 20)], [(7, 56, 24), (1, 45, 18)], [(10, 31, 26), (1, 27, 26)], [(10, 17, 26), (3, 20, 28), (1, 21, 28)]],
    [[(3, 139, 26), (1, 137, 28)], [(6, 66, 28), (1, 66, 30)], [(9, 36, 30), (1, 34, 32)], [(13, 19, 28), (1, 17, 32)]],
    [[(6, 84, 16), (1, 82, 16)], [(6, 70, 30), (1, 68, 30)], [(7, 35, 30), (3, 33, 28), (1, 32, 28)], [(13, 20, 30), (1, 20, 28)]],
    [[(5, 105, 20), (1, 94, 18)], [(6, 74, 32), (1, 71, 30)], [(11, 33, 28), (1, 34, 32)], [(13, 19, 28), (3, 16, 26)]],
    [[(4, 127, 24), (1, 126, 24)], [(7, 66, 28), (1, 66, 30)], [(12, 30, 24), (1, 24, 28), (1, 24, 30)], [(15, 19, 28), (1, 17, 32)]],
    [[(7, 84, 16), (1, 78, 16)], [(7, 70, 30), (1, 66, 28)], [(12, 33, 28), (1, 32, 30)], [(14, 21, 32), (1, 24, 28)]],
    [[(5, 117, 22), (1, 117, 24)], [(8, 66, 28), (1, 58, 26)], [(11, 38, 32), (1, 34, 32)], [(15, 20, 30), (2, 17, 26)]],
    [[(4, 148, 28), (1, 146, 28)], [(8, 68, 30), (1, 70, 24)], [(10, 36, 32), (3, 38, 28)], [(16, 19, 28), (3, 16, 26)]],
    [[(4, 126, 24), (2, 135, 26)], [(8, 70, 28), (2, 43, 26)], [(13, 32, 28), (2, 41, 30)], [(17, 19, 28), (3, 15, 26)]],
    [[(5, 136, 26), (1, 132, 24)], [(5, 67, 30), (4, 68, 28), (1, 69, 28)], [(14, 35, 30), (1, 32, 24)], [(18, 18, 26), (3, 16, 28), (1, 14, 28)]],
    [[(3, 142, 26), (3, 141, 28)], [(8, 70, 30), (1, 73, 32), (1, 74, 32)], [(12, 34, 30), (3, 34, 26), (1, 35, 28)], [(18, 21, 32), (1, 27, 30)]],
    [[(5, 116, 22), (2, 103, 20), (1, 102, 20)], [(9, 74, 32), (1, 74, 30)], [(14, 34, 28), (2, 32, 32), (1, 32, 30)], [(19, 21, 32), (1, 25, 26)]],
    [[(7, 116, 22), (1, 117, 22)], [(11, 65, 28), (1, 58, 24)], [(15, 38, 32), (1, 27, 28)], [(20, 20, 30), (1, 20, 32), (1, 21, 32)]],
    [[(6, 136, 26), (1, 130, 24)], [(11, 66, 28), (1, 62, 30)], [(14, 34, 28), (3, 34, 32), (1, 30, 30)], [(18, 20, 30), (3, 20, 28), (2, 15, 26)]],
    [[(5, 105, 20), (2, 115, 22), (2, 116, 22)], [(10, 75, 32), (1, 73, 32)], [(16, 38, 32), (1, 27, 28)], [(22, 19, 28), (2, 16, 30), (1, 19, 30)]],
    [[(6, 147, 28), (1, 146, 28)], [(11, 66, 28), (2, 65, 30)], [(18, 33, 28), (2, 33, 30)], [(22, 21, 32), (1, 28, 30)]],
    [[(6, 116, 22), (3, 125, 24)], [(11, 75, 32), (1, 68, 30)], [(13, 35, 28), (6, 34, 32), (1, 30, 30)], [(23, 21, 32), (1, 26, 30)]],
    [[(7, 105, 20), (4, 95, 18)], [(12, 67, 28), (1, 63, 30), (1, 62, 32)], [(21, 31, 26), (2, 33, 32)], [(23, 21, 32), (2, 24, 30)]],
    [[(10, 116, 22)], [(12, 74, 32), (1, 78, 30)], [(18, 37, 32), (1, 39, 30), (1, 41, 28)], [(25, 21, 32), (1, 27, 28)]],
    [[(5, 126, 24), (4, 115, 22), (1, 114, 22)], [(12, 67, 28), (2, 66, 32), (1, 68, 30)], [(21, 35, 30), (1, 39, 30)], [(26, 21, 32), (1, 28, 28)]],
    [[(9, 126, 24), (1, 117, 22)], [(13, 75, 32), (1, 68, 30)], [(20, 35, 30), (3, 35, 28)], [(27, 21, 32), (1, 28, 30)]],
    [[(9, 126, 24), (1, 137, 26)], [(13, 71, 30), (2, 68, 32)], [(20, 37, 32), (1, 39, 28), (1, 38, 28)], [(24, 20, 32), (5, 25, 28)]],
    [[(8, 147, 28), (1, 141, 28)], [(10, 73, 32), (4, 74, 30), (1, 73, 30)], [(16, 36, 32), (6, 39, 30), (1, 37, 30)], [(27, 21, 32), (3, 20, 26)]],
    [[(9, 137, 26), (1, 135, 26)], [(12, 70, 30), (4, 75, 32)], [(24, 35, 30), (1, 40, 28)], [(23, 20, 32), (8, 24, 30)]],
    [[(14, 95, 18), (1, 86, 18)], [(13, 73, 32), (3, 77, 30)], [(24, 35, 30), (2, 35, 28)], [(26, 21, 32), (5, 21, 30), (1, 23, 30)]],
    [[(9, 147, 28), (1, 142, 28)], [(10, 73, 30), (6, 70, 32), (1, 71, 32)], [(25, 35, 30), (2, 34, 26)], [(29, 21, 32), (4, 22, 30)]],
    [[(11, 126, 24), (1, 131, 24)], [(16, 74, 32), (1, 79, 30)], [(25, 38, 32), (1, 25, 30)], [(33, 21, 32), (1, 28, 28)]],
    [[(14, 105, 20), (1, 99, 18)], [(19, 65, 28), (1, 72, 28)], [(24, 37, 32), (2, 40, 30), (1, 41, 30)], [(31, 21, 32), (4, 24, 32)]],
    [[(10, 147, 28), (1, 151, 28)], [(15, 71, 30), (3, 71, 32), (1, 73, 32)], [(24, 37, 32), (3, 38, 30), (1, 39, 30)], [(36, 19, 30), (3, 29, 26)]],
    [[(15, 105, 20), (1, 99, 18)], [(19, 70, 30), (1, 64, 28)], [(27, 38, 32), (2, 25, 26)], [(38, 20, 30), (2, 18, 28)]],
    [[(14, 105, 20), (1, 113, 22), (1, 114, 22)], [(17, 67, 30), (3, 92, 32)], [(30, 35, 30), (1, 41, 30)], [(36, 21, 32), (1, 26, 30), (1, 27, 30)]],
    [[(11, 146, 28), (1, 146, 26)], [(20, 70, 30), (1, 60, 26)], [(29, 38, 32), (1, 24, 32)], [(40, 20, 30), (2, 17, 26)]],
    [[(3, 137, 26), (1, 136, 26), (10, 126, 24)], [(22, 65, 28), (1, 75, 30)], [(30, 37, 32), (1, 51, 30)], [(42, 20, 30), (1, 21, 30)]],
    [[(12, 126, 24), (2, 118, 22), (1, 116, 22)], [(19, 74, 32), (1, 74, 30), (1, 72, 28)], [(30, 38, 32), (2, 29, 30)], [(39, 20, 32), (2, 37, 26), (1, 38, 26)]],
    [[(12, 126, 24), (3, 136, 26)], [(21, 70, 30), (2, 65, 28)], [(34, 35, 30), (1, 44, 32)], [(42, 20, 30), (2, 19, 28), (2, 18, 28)]],
    [[(12, 126, 24), (3, 117, 22), (1, 116, 22)], [(25, 61, 26), (2, 62, 28)], [(34, 35, 30), (1, 40, 32), (1, 41, 32)], [(45, 20, 30), (1, 20, 32), (1, 21, 32)]],
    [[(15, 105, 20), (2, 115, 22), (2, 116, 22)], [(25, 65, 28), (1, 72, 28)], [(18, 35, 30), (17, 37, 32), (1, 50, 32)], [(42, 20, 30), (6, 19, 28), (1, 15, 28)]],
    [[(19, 105, 20), (1, 101, 20)], [(33, 51, 22), (1, 65, 22)], [(40, 33, 28), (1, 28, 28)], [(49, 20, 30), (1, 18, 28)]],
    [[(18, 105, 20), (2, 117, 22)], [(26, 65, 28), (1, 80, 30)], [(35, 35, 30), (3, 35, 28), (1, 36, 28)], [(52, 18, 28), (2, 38, 30)]],
    [[(26, 84, 16)], [(26, 70, 30)], [(45, 31, 26), (1, 9, 26)], [(52, 20, 30)]],
    [[(16, 126, 24), (1, 114, 22), (1, 115, 22)], [(23, 70, 30), (3, 65, 28), (1, 66, 28)], [(40, 35, 30), (1, 43, 30)], [(46, 20, 30), (7, 19, 28), (1, 16, 28)]],
    [[(19, 116, 22), (1, 105, 22)], [(20, 70, 30), (7, 66, 28), (1, 63, 28)], [(40, 35, 30), (1, 42, 32), (1, 43, 32)], [(54, 20, 30), (1, 19, 30)]],
    [[(17, 126, 24), (2, 115, 22)], [(24, 70, 30), (4, 74, 32)], [(48, 31, 26), (2, 18, 26)], [(54, 19, 28), (6, 15, 26), (1, 14, 26)]],
    [[(29, 84, 16)], [(29, 70, 30)], [(6, 34, 30), (3, 36, 30), (38, 33, 28)], [(58, 20, 30)]],
    [[(16, 147, 28), (1, 149, 28)], [(31, 66, 28), (1, 37, 26)], [(48, 33, 28), (1, 23, 26)], [(53, 20, 30), (6, 19, 28), (1, 17, 28)]],
    [[(20, 115, 22), (2, 134, 24)], [(29, 66, 28), (2, 56, 26), (2, 57, 26)], [(45, 36, 30), (2, 15, 28)], [(59, 20, 30), (2, 21, 32)]],
    [[(17, 147, 28), (1, 134, 26)], [(26, 70, 30), (5, 75, 32)], [(47, 35, 30), (1, 48, 32)], [(64, 18, 28), (2, 33, 30), (1, 35, 30)]],
    [[(22, 115, 22), (1, 133, 24)], [(33, 65, 28), (1, 74, 28)], [(43, 36, 30), (5, 27, 28), (1, 30, 28)], [(57, 20, 30), (5, 21, 32), (1, 24, 32)]],
    [[(18, 136, 26), (2, 142, 26)], [(33, 66, 28), (2, 49, 26)], [(48, 35, 30), (2, 38, 28)], [(64, 20, 30), (1, 20, 32)]],
    [[(19, 126, 24), (2, 135, 26), (1, 136, 26)], [(32, 66, 28), (2, 55, 26), (2, 56, 26)], [(49, 36, 30), (2, 18, 32)], [(65, 18, 28), (5, 27, 30), (1, 29, 30)]],
    [[(20, 137, 26), (1, 130, 26)], [(30, 75, 32), (2, 71, 32)], [(46, 35, 30), (6, 39, 32)], [(3, 12, 30), (70, 19, 28)]],
    [[(20, 147, 28)], [(35, 70, 30)], [(49, 35, 30), (5, 35, 28)], [(70, 20, 30)]],
    [[(21, 136, 26), (1, 155, 28)], [(34, 70, 30), (1, 64, 28), (1, 65, 28)], [(54, 35, 30), (1, 45, 30)], [(68, 20, 30), (3, 18, 28), (1, 19, 28)]],
    [[(19, 126, 24), (5, 115, 22), (1, 114, 22)], [(33, 70, 30), (3, 65, 28), (1, 64, 28)], [(52, 35, 30), (3, 41, 32), (1, 40, 32)], [(67, 20, 30), (5, 21, 32), (1, 24, 32)]],
    [[(2, 150, 28), (21, 136, 26)], [(32, 70, 30), (6, 65, 28)], [(52, 38, 32), (2, 27, 32)], [(73, 20, 30), (2, 22, 32)]],
    [[(21, 126, 24), (4, 136, 26)], [(30, 74, 32), (6, 73, 30)], [(54, 35, 30), (4, 40, 32)], [(75, 20, 30), (1, 20, 28)]],
    [[(30, 105, 20), (1, 114, 22)], [(3, 45, 22), (55, 47, 20)], [(2, 26, 26), (62, 33, 28)], [(79, 18, 28), (4, 33, 30)]],
]

# GB/T 21049-2022 表 A.1 校正图形参数 (m x m 区域划分, r/k 为区域尺寸)
ALIGN_R = [0, 0, 0, 15, 15, 17, 18, 19, 20, 21, 15, 15, 15, 17, 17, 19, 19, 19, 19, 21, 21, 17, 16, 18, 17, 19, 18, 20, 19, 21, 20, 17, 19, 17, 19, 17, 19, 21, 19, 21, 18, 20, 17, 19, 21, 18, 20, 22, 17, 19, 15, 17, 19, 21, 17, 19, 21, 18, 20, 15, 17, 19, 21, 16, 18, 17, 19, 21, 15, 17, 19, 21, 15, 17, 18, 20, 22, 15, 17, 19, 21, 23, 17, 19]
ALIGN_K = [0, 0, 0, 14, 16, 16, 17, 18, 19, 20, 14, 15, 16, 16, 17, 17, 18, 19, 20, 20, 21, 16, 17, 17, 18, 18, 19, 19, 20, 20, 21, 17, 17, 18, 18, 19, 19, 19, 20, 20, 17, 17, 18, 18, 18, 19, 19, 19, 17, 17, 18, 18, 18, 18, 19, 19, 19, 17, 17, 18, 18, 18, 18, 19, 19, 17, 17, 17, 18, 18, 18, 18, 19, 19, 17, 17, 17, 18, 18, 18, 18, 18, 17, 17]
ALIGN_M = [0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10]

# GB/T 21049-2022 表 18 URI 模式字符集 (生成端未使用, 仅解码端需要)
URI_SET_A = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", ".", "/", "-", "_", "~", ":", "@", "?", "#", "=", "+", "$", "&", "http://", "https://", "ftp://", "mailto:", "ldap://", "tel:", "urn:", "www.", ".com", ".net", ".gov", ".org", ".cn"]
URI_SET_B = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "!", "*", "(", ")", ",", "{", "}", "|", "\\", "^", "[", "]", "`", "<", ">", "%", "\"", ";", ".htm", ".html", ".asp", ".aspx", ".php", ".jsp", "gtin", "ser", "bat", "exp", "search", "id", ".jp", ".it", ".de", ".br", ".fr", "gs1"]
URI_SET_C = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z", "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "$", "-", "_", ".", "+", "!", "*", "(", ")", ",", "{", "}", "|", "\\", "^", "~", "[", "]", "`", "<", ">", "#", "%", "\"", ";", "/", "?", ":", "@", "&", "=", "http://", "https://", "ftp://", "mailto:", "ldap://", "tel:", "urn:", "www.", ".com", ".net", ".gov", ".org", ".cn", ".htm", ".html", ".asp", ".aspx", ".php", ".jsp", "gtin", "ser", "bat", "exp", "search", "id", ".jp", ".it", ".de", ".br", ".fr", "gs1", "search"]


@lru_cache(maxsize=None)
def rs_structure(version: int, level: int) -> dict:
    """返回某版本/纠错等级下的码块结构。"""
    batches = RS_BLOCKS[version - 1][level]
    blocks = [tuple(batch) for batch in batches]
    total = sum(count * (data + corr) for count, data, corr in blocks)
    data_codewords = sum(count * data for count, data, _ in blocks)
    correction_codewords = sum(count * corr for count, _, corr in blocks)
    return {
        "version": version,
        "level": level,
        "batches": blocks,
        "block_count": sum(count for count, _, _ in blocks),
        "total_codewords": total,
        "data_codewords": data_codewords,
        "correction_codewords": correction_codewords,
        "data_bits": data_codewords * 8,
    }


def expand_rs_blocks(version: int, level: int):
    blocks = []
    for batch_index, (count, data, corr) in enumerate(RS_BLOCKS[version - 1][level]):
        for _ in range(count):
            blocks.append((batch_index, data, corr))
    return blocks


# --------------------------------------------------------------------------
# 符号模板: 四角定位图形 / 分隔符 / 校正图形 / 辅助校正图形 / 功能信息区
# --------------------------------------------------------------------------
FINDER_PATTERN = (
    (1, 1, 1, 1, 1, 1, 1),
    (1, 0, 0, 0, 0, 0, 0),
    (1, 0, 1, 1, 1, 1, 1),
    (1, 0, 1, 0, 0, 0, 0),
    (1, 0, 1, 0, 1, 1, 1),
    (1, 0, 1, 0, 1, 1, 1),
    (1, 0, 1, 0, 1, 1, 1),
)
FINDER_SIZE = 7


def _finder_value(rotation: int, x: int, y: int) -> int:
    if rotation == 1:
        source_x, source_y = y, FINDER_SIZE - 1 - x
    elif rotation == 2:
        source_x, source_y = FINDER_SIZE - 1 - x, FINDER_SIZE - 1 - y
    else:
        source_x, source_y = x, y
    return FINDER_PATTERN[source_y][source_x]


def alignment_spans(version: int):
    """把符号切成 m x m 区域的边界尺寸(GB/T 21049-2022 表 A.1)。"""
    m = ALIGN_M[version - 1]
    if m is None:
        return []
    r = ALIGN_R[version - 1]
    k = ALIGN_K[version - 1]
    repeated = k if k is not None else dimension_for_version(version) - r
    return [repeated] * m + [r - 1]


class Template:
    __slots__ = ("version", "dimension", "modules", "functions")

    def __init__(self, version, dimension, modules, functions):
        self.version = version
        self.dimension = dimension
        self.modules = modules
        self.functions = functions


@lru_cache(maxsize=None)
def build_template(version: int) -> Template:
    dimension = dimension_for_version(version)
    modules = [0] * (dimension * dimension)
    functions = [0] * (dimension * dimension)

    def plot(column, row, value):
        if 0 <= column < dimension and 0 <= row < dimension:
            if not functions[row * dimension + column]:
                modules[row * dimension + column] = 1 if value else 0
                functions[row * dimension + column] = 1

    def reserve(column, row):
        functions[row * dimension + column] = 1

    def place_finder(origin_column, origin_row, rotation):
        for y in range(FINDER_SIZE):
            for x in range(FINDER_SIZE):
                column, row = origin_column + x, origin_row + y
                modules[row * dimension + column] = _finder_value(rotation, x, y)
                functions[row * dimension + column] = 1

    far = dimension - FINDER_SIZE
    place_finder(0, 0, 0)
    place_finder(far, 0, 1)
    place_finder(0, far, 1)
    place_finder(far, far, 2)

    separator_far = dimension - 8
    for offset in range(8):
        reserve(offset, 7)
        reserve(7, offset)
        reserve(dimension - 1 - offset, 7)
        reserve(separator_far, offset)
        reserve(offset, separator_far)
        reserve(7, dimension - 1 - offset)
        reserve(dimension - 1 - offset, separator_far)
        reserve(separator_far, dimension - 1 - offset)

    primary, secondary = function_information_coordinates(dimension)
    for column, row in primary:
        reserve(column, row)
    for column, row in secondary:
        reserve(column, row)

    spans = alignment_spans(version)
    if spans:
        def place_assistant_pattern(center_column, center_row):
            for row_offset in (-1, 0, 1):
                for column_offset in (-1, 0, 1):
                    plot(
                        center_column + column_offset,
                        center_row + row_offset,
                        column_offset == 0 and row_offset == 0,
                    )

        def place_alignment_pattern(right, top, width, height):
            plot(right, top, True)
            plot(right - 1, top + 1, False)
            for offset in range(1, width + 1):
                plot(right - offset, top, True)
                plot(right - offset - 1, top + 1, False)
            for offset in range(1, height):
                plot(right, top + offset, True)
                plot(right - 1, top + offset + 1, False)

        m = len(spans) - 1
        boundary_count = len(spans) + 1

        row = 0
        for region in range(boundary_count):
            if (region % 2 == 0 and m % 2 == 1) or (region % 2 == 1 and m % 2 == 0):
                place_assistant_pattern(0, row)
            if region % 2 == 1:
                place_assistant_pattern(dimension - 1, row)
            if region < len(spans):
                row += spans[region]

        column = dimension - 1
        for region in range(boundary_count):
            if (region % 2 == 0 and m % 2 == 1) or (region % 2 == 1 and m % 2 == 0):
                place_assistant_pattern(column, dimension - 1)
            if region % 2 == 1:
                place_assistant_pattern(column, 0)
            if region < len(spans):
                column -= spans[region]

        top = 0
        for row_region in range(boundary_count):
            right = dimension - 1
            should_plot = row_region % 2 == 0
            height = spans[min(row_region, len(spans) - 1)]
            for column_region in range(boundary_count):
                width = spans[min(column_region, len(spans) - 1)]
                if should_plot and not (top == 0 and right == dimension - 1):
                    place_alignment_pattern(right, top, width, height)
                should_plot = not should_plot
                if column_region < len(spans):
                    right -= spans[column_region]
            if row_region < len(spans):
                top += spans[row_region]

    return Template(version, dimension, modules, functions)


def template_modules(version: int):
    template = build_template(version)
    return list(template.modules), list(template.functions)


# --------------------------------------------------------------------------
# 功能信息: 版本 + 纠错等级 + 掩模号, GF(16) 3 数据位 + 4 纠错位, 两份副本
# --------------------------------------------------------------------------
FUNCTION_INFO_DATA_SYMBOLS = 3
FUNCTION_INFO_CORRECTION_SYMBOLS = 4
FUNCTION_INFO_FILLER_BITS = 6


def function_information_coordinates(dimension: int):
    version_for_dimension(dimension)
    near, far = 8, dimension - 9
    primary, secondary = [], []
    for index in range(9):
        primary.append((index, near))
        secondary.append((dimension - 1 - index, far))
    for index in range(1, 9):
        primary.append((near, near - index))
        secondary.append((far, far + index))
    for index in range(9):
        primary.append((far, index))
        secondary.append((near, dimension - 1 - index))
    for index in range(1, 9):
        primary.append((far + index, near))
        secondary.append((near - index, far))
    return primary, secondary


def encode_function_information(version: int, level: int, mask: int):
    value = ((version + 20) << 4) | (level << 2) | mask
    data = [(value >> 8) & 0x0F, (value >> 4) & 0x0F, value & 0x0F]
    return rs_encode(GF16, data, FUNCTION_INFO_CORRECTION_SYMBOLS, 1)


def decode_function_information(codewords):
    corrected, errors = rs_decode(GF16, codewords, FUNCTION_INFO_CORRECTION_SYMBOLS, 1)
    high, low, settings = corrected[:FUNCTION_INFO_DATA_SYMBOLS]
    version = ((high << 4) | low) - 20
    if not MIN_VERSION <= version <= MAX_VERSION:
        raise FunctionInformationError(f"功能信息里的版本号 {version} 非法")
    return {
        "version": version,
        "level": settings >> 2,
        "mask": settings & 0b11,
        "corrected": errors,
    }


def _bits_to_symbols(bits):
    symbols = [0] * (FUNCTION_INFO_DATA_SYMBOLS + FUNCTION_INFO_CORRECTION_SYMBOLS)
    for index, bit in enumerate(bits):
        symbols[index >> 2] = (symbols[index >> 2] << 1) | bit
    return symbols


def read_function_information(modules, dimension):
    primary, secondary = function_information_coordinates(dimension)
    copies = []
    for coordinates in (primary, secondary):
        bits = [modules[row * dimension + column] for column, row in coordinates]
        copies.append(bits)
    results, errors = [], []
    for bits in copies:
        try:
            results.append(decode_function_information(_bits_to_symbols(bits[:28])))
        except HanXinError as error:
            errors.append(error)
            results.append(None)
    if results[0] and results[1]:
        if (results[0]["version"], results[0]["level"], results[0]["mask"]) != (
            results[1]["version"], results[1]["level"], results[1]["mask"]
        ):
            raise FunctionInformationError("两份功能信息副本解出的参数不一致")
    selected = results[0] or results[1]
    if selected is None:
        raise FunctionInformationError("两份功能信息副本都无法解码") from errors[0]
    selected = dict(selected)
    selected["sources"] = "both" if results[0] and results[1] else ("primary" if results[0] else "secondary")
    return selected


# --------------------------------------------------------------------------
# 掩模 (GB/T 21049-2022 表 14), i/j 为 1 起的行/列号
# --------------------------------------------------------------------------
def mask_bit(mask: int, row: int, column: int) -> bool:
    i, j = row + 1, column + 1
    if mask == 0:
        return False
    if mask == 1:
        return (i + j) % 2 == 0
    if mask == 2:
        return (((i + j) % 3) + (j % 3)) % 2 == 0
    if mask == 3:
        return ((i % j) + (j % i) + (i % 3) + (j % 3)) % 2 == 0
    raise ValueError("掩模号必须是 0..3")


# --------------------------------------------------------------------------
# 数据区排布: 按行优先填入非功能模块, 码字按 13 路栅栏交织
# --------------------------------------------------------------------------
PICKET_FENCE_STRIDE = 13


@lru_cache(maxsize=None)
def data_placement(version: int):
    template = build_template(version)
    total_codewords = rs_structure(version, 0)["total_codewords"]
    bit_count = total_codewords * 8
    offsets = [index for index, flag in enumerate(template.functions) if not flag]
    if len(offsets) < bit_count:
        raise HanXinError(f"版本 {version} 的数据区模块不足")
    return tuple(offsets[:bit_count]), tuple(offsets[bit_count:])


def to_picket_fence(sequential):
    ordered = [0] * len(sequential)
    position = 0
    for start in range(PICKET_FENCE_STRIDE):
        for index in range(start, len(sequential), PICKET_FENCE_STRIDE):
            ordered[position] = sequential[index]
            position += 1
    return ordered


def from_picket_fence(placed):
    sequential = [0] * len(placed)
    position = 0
    for start in range(PICKET_FENCE_STRIDE):
        for index in range(start, len(placed), PICKET_FENCE_STRIDE):
            sequential[index] = placed[position]
            position += 1
    return sequential


def correct_blocks(placed_codewords, version: int, level: int):
    sequential = from_picket_fence(placed_codewords)
    structure = rs_structure(version, level)
    if len(sequential) != structure["total_codewords"]:
        raise HanXinError("栅栏交织后的码字数与版本容量不符")
    data_codewords = bytearray()
    reports = []
    offset = 0
    corrected_errors = 0
    for batch_index, data_count, correction_count in expand_rs_blocks(version, level):
        block = sequential[offset:offset + data_count + correction_count]
        offset += data_count + correction_count
        corrected, errors = rs_decode(GF256, block, correction_count, 1)
        corrected_errors += errors
        data_codewords.extend(corrected[:data_count])
        reports.append({
            "batch": batch_index,
            "data_codewords": data_count,
            "correction_codewords": correction_count,
            "corrected_errors": errors,
        })
    return bytes(data_codewords), corrected_errors, reports


# --------------------------------------------------------------------------
# 位流读写
# --------------------------------------------------------------------------
class BitReader:
    __slots__ = ("data", "bit_length", "position")

    def __init__(self, data, bit_length=None):
        self.data = bytes(data)
        self.bit_length = len(self.data) * 8 if bit_length is None else bit_length
        self.position = 0

    @property
    def available(self) -> int:
        return self.bit_length - self.position

    def read(self, count: int) -> int:
        if count > self.available:
            raise InvalidBitStreamError(f"需要 {count} 位, 仅剩 {self.available} 位")
        value = 0
        for _ in range(count):
            byte = self.data[self.position >> 3]
            value = (value << 1) | ((byte >> (7 - (self.position & 7))) & 1)
            self.position += 1
        return value

    def peek(self, count: int) -> int:
        position = self.position
        try:
            return self.read(count)
        finally:
            self.position = position

    def skip(self, count: int) -> None:
        if count > self.available:
            raise InvalidBitStreamError("跳过位数超出可用位数")
        self.position += count

    def has_only_zeroes(self) -> bool:
        position = self.position
        try:
            while self.available > 0:
                if self.read(min(32, self.available)):
                    return False
            return True
        finally:
            self.position = position


class BitWriter:
    __slots__ = ("bits",)

    def __init__(self):
        self.bits = []

    def write(self, value: int, count: int) -> None:
        if value < 0 or value >> count:
            raise ValueError(f"数值 {value} 无法放入 {count} 位")
        for shift in range(count - 1, -1, -1):
            self.bits.append((value >> shift) & 1)

    def __len__(self) -> int:
        return len(self.bits)


# --------------------------------------------------------------------------
# 信息位流: 模式解析 (GB/T 21049-2022 第 6 章)
# --------------------------------------------------------------------------
MODE_NUMERIC = "numeric"
MODE_TEXT = "text"
MODE_BINARY = "binary"
MODE_CHINESE_ONE = "chinese-region-one"
MODE_CHINESE_TWO = "chinese-region-two"
MODE_GB18030_TWO = "gb18030-two-byte"
MODE_GB18030_FOUR = "gb18030-four-byte"
MODE_ECI = "eci"
MODE_UNICODE = "unicode"
MODE_GS1 = "gs1"
MODE_URI = "uri"

SHORT_MODES = {
    0b0001: MODE_NUMERIC,
    0b0010: MODE_TEXT,
    0b0011: MODE_BINARY,
    0b0100: MODE_CHINESE_ONE,
    0b0101: MODE_CHINESE_TWO,
    0b0110: MODE_GB18030_TWO,
    0b0111: MODE_GB18030_FOUR,
    0b1000: MODE_ECI,
    0b1001: MODE_UNICODE,
}
LONG_MODES = {0b11100001: MODE_GS1, 0b11100010: MODE_URI}

NUMERIC_TERMINATORS = {1021: 1, 1022: 2, 1023: 3}
TEXT_SWITCH = 62
TEXT_TERMINATOR = 63
CHINESE_SWITCH = 0xFFE
CHINESE_TERMINATOR = 0xFFF
GB18030_TWO_BYTE_TERMINATOR = 0x7FFF
GS1_FNC1 = 1000


def read_mode_indicator(reader: BitReader) -> str:
    if reader.available < 4:
        raise InvalidBitStreamError("模式指示符被截断")
    prefix = reader.read(4)
    if prefix in SHORT_MODES:
        return SHORT_MODES[prefix]
    if prefix != 0b1110:
        raise InvalidBitStreamError(f"保留的模式指示符 {prefix:04b}")
    if reader.available < 4:
        raise InvalidBitStreamError("8 位模式指示符被截断")
    indicator = (prefix << 4) | reader.read(4)
    if indicator not in LONG_MODES:
        raise InvalidBitStreamError(f"保留的模式指示符 {indicator:08b}")
    return LONG_MODES[indicator]


def _ascii_bytes(text: str) -> bytes:
    return bytes(ord(character) & 0xFF for character in text)


def _decode_gb18030(data):
    try:
        return bytes(data).decode("gb18030")
    except UnicodeDecodeError:
        return None


def _decode_text_fallback(data):
    """字节段没有声明字符集时, 依次尝试 UTF-8 / GB 18030 还原文本。"""
    for encoding in ("utf-8", "gb18030"):
        try:
            return bytes(data).decode(encoding)
        except UnicodeDecodeError:
            continue
    return None


def _text_one_byte(value: int) -> int:
    if value <= 9:
        return 0x30 + value
    if value <= 35:
        return 0x41 + value - 10
    return 0x61 + value - 36


def _text_two_byte(value: int) -> int:
    if value <= 27:
        return value
    if value <= 43:
        return 0x20 + value - 28
    if value <= 50:
        return 0x3A + value - 44
    if value <= 56:
        return 0x5B + value - 51
    return 0x7B + value - 57


def read_numeric_segment(reader: BitReader) -> dict:
    groups = []
    while True:
        if reader.available < 10:
            raise InvalidBitStreamError("数字段缺少完整的终止符")
        value = reader.read(10)
        final_digits = NUMERIC_TERMINATORS.get(value)
        if final_digits is not None:
            if not groups:
                raise InvalidBitStreamError("数字段终止符出现在任何数字之前")
            text = "".join(
                str(group).zfill(final_digits if index == len(groups) - 1 else 3)
                for index, group in enumerate(groups)
            )
            return {
                "mode": MODE_NUMERIC,
                "text": text,
                "bytes": _ascii_bytes(text),
                "character_count": len(text),
            }
        if value > 999:
            raise InvalidBitStreamError(f"保留的数字组 {value}")
        groups.append(value)


def read_text_segment(reader: BitReader) -> dict:
    values = []
    submode = 1
    needs_character = False
    while True:
        if reader.available < 6:
            raise InvalidBitStreamError("文本段缺少完整的终止符")
        value = reader.read(6)
        if value == TEXT_TERMINATOR:
            if not values:
                raise InvalidBitStreamError("文本段终止符出现在任何字符之前")
            if needs_character:
                raise InvalidBitStreamError("文本子模式切换符后没有字符")
            break
        if value == TEXT_SWITCH:
            if needs_character:
                raise InvalidBitStreamError("连续的文本子模式切换符")
            submode = 2 if submode == 1 else 1
            needs_character = True
            continue
        values.append(_text_one_byte(value) if submode == 1 else _text_two_byte(value))
        needs_character = False
    data = bytes(values)
    return {
        "mode": MODE_TEXT,
        "text": data.decode("latin-1"),
        "bytes": data,
        "character_count": len(data),
    }


def read_binary_segment(reader: BitReader) -> dict:
    if reader.available < 13:
        raise InvalidBitStreamError("字节段缺少字节数")
    count = reader.read(13)
    if count == 0:
        raise InvalidBitStreamError("字节段长度必须为正")
    if reader.available < count * 8:
        raise InvalidBitStreamError("字节段声明的字节数超出剩余位流")
    data = bytes(reader.read(8) for _ in range(count))
    return {
        "mode": MODE_BINARY,
        "bytes": data,
        "text": _decode_text_fallback(data),
        "character_count": count,
    }


def _chinese_region_one_bytes(value: int):
    if value < 3760:
        return bytes(((value // 94) + 0xB0, (value % 94) + 0xA1))
    if value < 4042:
        offset = value - 3760
        return bytes(((offset // 94) + 0xA1, (offset % 94) + 0xA1))
    if value < 4074:
        return bytes((0xA8, value - 4042 + 0xA1))
    raise InvalidBitStreamError(f"保留的常用汉字一区值 {value}")


def _chinese_region_two_bytes(value: int):
    if value >= 3008:
        raise InvalidBitStreamError(f"保留的常用汉字二区值 {value}")
    return bytes(((value // 94) + 0xD8, (value % 94) + 0xA1))


def read_common_chinese_segment(reader: BitReader, initial_region: int) -> dict:
    region = initial_region
    values = []
    output = bytearray()
    while True:
        if reader.available < 12:
            raise InvalidBitStreamError("常用汉字段缺少完整的终止符")
        value = reader.read(12)
        if value == CHINESE_TERMINATOR:
            break
        if value == CHINESE_SWITCH:
            region = 2 if region == 1 else 1
            continue
        if region == 1:
            output += _chinese_region_one_bytes(value)
        else:
            output += _chinese_region_two_bytes(value)
        values.append((region, value))
    if not values:
        raise InvalidBitStreamError("常用汉字段终止符出现在任何汉字之前")
    return {
        "mode": MODE_CHINESE_ONE if initial_region == 1 else MODE_CHINESE_TWO,
        "bytes": bytes(output),
        "text": _decode_gb18030(output),
        "character_count": len(values),
    }


def read_gb18030_two_byte_segment(reader: BitReader) -> dict:
    output = bytearray()
    count = 0
    while True:
        if reader.available < 15:
            raise InvalidBitStreamError("GB 18030 双字节段缺少完整的终止符")
        value = reader.read(15)
        if value == GB18030_TWO_BYTE_TERMINATOR:
            break
        if value >= 23940:
            raise InvalidBitStreamError(f"保留的 GB 18030 双字节值 {value}")
        first = value // 190 + 0x81
        trail = value % 190
        second = trail + 0x40 if trail < 63 else trail + 0x41
        output += bytes((first, second))
        count += 1
    if count == 0:
        raise InvalidBitStreamError("GB 18030 双字节段终止符出现在任何字符之前")
    return {
        "mode": MODE_GB18030_TWO,
        "bytes": bytes(output),
        "text": _decode_gb18030(output),
        "character_count": count,
    }


def read_gb18030_four_byte_segment(reader: BitReader) -> dict:
    if reader.available < 21:
        raise InvalidBitStreamError("GB 18030 四字节段被截断")
    value = reader.read(21)
    if value >= 1_587_600:
        raise InvalidBitStreamError(f"保留的 GB 18030 四字节值 {value}")
    first, remainder = divmod(value, 12_600)
    second, remainder = divmod(remainder, 1_260)
    third, fourth = divmod(remainder, 10)
    data = bytes((first + 0x81, second + 0x30, third + 0x81, fourth + 0x30))
    return {
        "mode": MODE_GB18030_FOUR,
        "bytes": data,
        "text": _decode_gb18030(data),
        "character_count": 1,
    }


def read_eci_segment(reader: BitReader) -> dict:
    if reader.available < 8:
        raise InvalidBitStreamError("ECI 指派号被截断")
    prefix = reader.peek(2)
    if (prefix & 0b10) == 0:
        assignment = reader.read(8) & 0x7F
    elif prefix == 0b10:
        if reader.available < 16:
            raise InvalidBitStreamError("16 位 ECI 指派号被截断")
        assignment = reader.read(16) & 0x3FFF
    else:
        if reader.available < 24:
            raise InvalidBitStreamError("24 位 ECI 指派号被截断")
        assignment = reader.read(24) & 0x1F_FFFF
        if assignment > 999_999:
            raise InvalidBitStreamError(f"保留的 ECI 指派号 {assignment}")
    return {
        "mode": MODE_ECI,
        "bytes": b"",
        "text": None,
        "character_count": 0,
        "assignment_number": assignment,
    }


def _read_unicode_group_count(reader: BitReader) -> int:
    if reader.available < 4:
        raise InvalidBitStreamError("Unicode 字节模式计数被截断")
    prefix = reader.peek(4)
    if (prefix & 0b1000) == 0:
        bit_length, payload_bits = 4, 3
    elif (prefix & 0b1100) == 0b1000:
        bit_length, payload_bits = 8, 6
    elif (prefix & 0b1110) == 0b1100:
        bit_length, payload_bits = 12, 9
    elif prefix == 0b1110:
        bit_length, payload_bits = 16, 12
    else:
        bit_length, payload_bits = 20, 15
    if reader.available < bit_length:
        raise InvalidBitStreamError(f"{bit_length} 位的 Unicode 计数被截断")
    value = reader.read(bit_length) & ((1 << payload_bits) - 1)
    if value == 0:
        raise InvalidBitStreamError("Unicode 字节模式计数必须为正")
    return value


def read_unicode_segment(reader: BitReader) -> dict:
    output = bytearray()
    groups = []
    while True:
        if reader.available < 4:
            raise InvalidBitStreamError("Unicode 段缺少完整的终止符")
        byte_mode = reader.read(4)
        if byte_mode == 0b1111:
            break
        if not 1 <= byte_mode <= 4:
            raise InvalidBitStreamError(f"保留的 Unicode 字节模式 {byte_mode}")
        count = _read_unicode_group_count(reader)
        if reader.available < byte_mode * 12:
            raise InvalidBitStreamError("Unicode 字节模式头被截断")
        widths = [reader.read(4) for _ in range(byte_mode)]
        if any(width > 8 for width in widths):
            raise InvalidBitStreamError("Unicode 差值位宽超过 8 位")
        minima = [reader.read(8) for _ in range(byte_mode)]
        for _ in range(count):
            values = []
            for byte_index in range(byte_mode):
                width = widths[byte_index]
                if reader.available < width:
                    raise InvalidBitStreamError("Unicode 差值数据被截断")
                value = minima[byte_index] + reader.read(width)
                if value > 0xFF:
                    raise InvalidBitStreamError("Unicode 差分字节超过 255")
                values.append(value)
            output += bytes(values)
        groups.append({"byte_mode": byte_mode, "count": count, "widths": widths, "minima": minima})
    if not groups:
        raise InvalidBitStreamError("Unicode 段终止符出现在任何字节组之前")
    try:
        text = bytes(output).decode("utf-8")
    except UnicodeDecodeError as error:
        raise InvalidBitStreamError("Unicode 段包含非法 UTF-8") from error
    return {
        "mode": MODE_UNICODE,
        "bytes": bytes(output),
        "text": text,
        "character_count": len(text),
    }


def _read_gs1_numeric_segment(reader: BitReader) -> dict:
    tokens = []
    final_digits = None
    while True:
        if reader.available < 10:
            raise InvalidBitStreamError("GS1 数字段缺少终止符")
        value = reader.read(10)
        digits = NUMERIC_TERMINATORS.get(value)
        if digits is not None:
            final_digits = digits
            break
        if value == GS1_FNC1:
            tokens.append(None)
        elif value <= 999:
            tokens.append(value)
        else:
            raise InvalidBitStreamError(f"保留的 GS1 数字值 {value}")
    numeric_tokens = [token for token in tokens if token is not None]
    if not numeric_tokens:
        raise InvalidBitStreamError("GS1 数字段终止符出现在任何数字之前")
    last = numeric_tokens[-1]
    text = ""
    for token in tokens:
        if token is None:
            text += "\x1d"
        else:
            text += str(token).zfill(final_digits if token is last else 3)
    return {
        "mode": MODE_NUMERIC,
        "text": text,
        "bytes": _ascii_bytes(text),
        "character_count": len(text),
        "gs1": True,
    }


def read_gs1_segment(reader: BitReader) -> dict:
    segments = []
    while True:
        if reader.available < 8:
            raise InvalidBitStreamError("GS1 段缺少终止符")
        if reader.peek(8) == 0xFF:
            reader.skip(8)
            break
        mode = read_mode_indicator(reader)
        if mode == MODE_NUMERIC:
            segment = _read_gs1_numeric_segment(reader)
        elif mode == MODE_TEXT:
            segment = dict(read_text_segment(reader), gs1=True)
        else:
            raise InvalidBitStreamError(f"GS1 内不允许 {mode} 模式")
        segments.append(segment)
    if not segments:
        raise InvalidBitStreamError("GS1 段终止符出现在任何数据段之前")
    data = bytearray()
    for segment in segments:
        data += segment["bytes"]
    return {
        "mode": MODE_GS1,
        "bytes": bytes(data),
        "text": "".join(segment["text"] or "" for segment in segments),
        "character_count": len(data),
        "segments": segments,
        "gs1": True,
    }


def _uri_definition(name: str):
    if name == "A":
        return URI_SET_A, 6, 62, 63
    if name == "B":
        return URI_SET_B, 6, 62, 63
    return URI_SET_C, 7, None, 127


def _read_uri_character_set(reader: BitReader, initial: str) -> str:
    name = initial
    text = ""
    while True:
        values, width, switch_value, end_value = _uri_definition(name)
        if reader.available < width:
            raise InvalidBitStreamError(f"URI-{name} 字符数据被截断")
        value = reader.read(width)
        if value == end_value:
            return text
        if name == "C" and value == 125:
            name = "A"
            continue
        if name == "C" and value == 126:
            name = "B"
            continue
        if value == switch_value:
            name = "B" if name == "A" else "A"
            continue
        if value >= len(values):
            raise InvalidBitStreamError(f"保留的 URI-{name} 值 {value}")
        text += values[value]


def read_uri_segment(reader: BitReader) -> dict:
    text = ""
    groups = []
    while True:
        if reader.available < 3:
            raise InvalidBitStreamError("URI 段缺少完整的终止符")
        indicator = reader.read(3)
        if indicator == 0b111:
            break
        if 1 <= indicator <= 3:
            name = chr(0x40 + indicator)
            value = _read_uri_character_set(reader, name)
            text += value
            groups.append({"set": name, "text": value})
            continue
        if indicator == 0b100:
            if reader.available < 6:
                raise InvalidBitStreamError("URI 百分号编码计数被截断")
            count = reader.read(6)
            if count == 0:
                raise InvalidBitStreamError("URI 百分号编码计数必须为正")
            if reader.available < count * 8:
                raise InvalidBitStreamError("URI 百分号编码字节被截断")
            escaped = ""
            for _ in range(count):
                escaped += f"%{reader.read(8):02X}"
            text += escaped
            groups.append({"set": "percent", "text": escaped})
            continue
        raise InvalidBitStreamError(f"保留的 URI 字符集指示符 {indicator}")
    if not groups:
        raise InvalidBitStreamError("URI 段终止符出现在任何字符之前")
    return {
        "mode": MODE_URI,
        "bytes": text.encode("utf-8"),
        "text": text,
        "character_count": len(text),
        "groups": groups,
        "uri": True,
    }


SEGMENT_READERS = {
    MODE_NUMERIC: read_numeric_segment,
    MODE_TEXT: read_text_segment,
    MODE_BINARY: read_binary_segment,
    MODE_CHINESE_ONE: lambda reader: read_common_chinese_segment(reader, 1),
    MODE_CHINESE_TWO: lambda reader: read_common_chinese_segment(reader, 2),
    MODE_GB18030_TWO: read_gb18030_two_byte_segment,
    MODE_GB18030_FOUR: read_gb18030_four_byte_segment,
    MODE_ECI: read_eci_segment,
    MODE_UNICODE: read_unicode_segment,
    MODE_GS1: read_gs1_segment,
    MODE_URI: read_uri_segment,
}


def read_payload(data) -> dict:
    reader = data if isinstance(data, BitReader) else BitReader(data)
    segments = []
    while reader.available > 0:
        if reader.has_only_zeroes():
            reader.skip(reader.available)
            break
        mode = read_mode_indicator(reader)
        segments.append(SEGMENT_READERS[mode](reader))

    data_segments = [segment for segment in segments if segment["mode"] != MODE_ECI]
    if not data_segments:
        raise InvalidBitStreamError("信息位流里没有数据段")
    payload = bytearray()
    for segment in data_segments:
        payload += segment["bytes"]
    has_text = all(segment["text"] is not None for segment in data_segments)
    result = {
        "bytes": bytes(payload),
        "text": "".join(segment["text"] for segment in data_segments) if has_text else None,
        "segments": segments,
        "bits_read": reader.position,
        "eci_used": any(segment["mode"] == MODE_ECI for segment in segments),
        "gs1": any(segment["mode"] == MODE_GS1 for segment in segments),
        "uri": any(segment["mode"] == MODE_URI for segment in segments),
        "unicode": any(segment["mode"] == MODE_UNICODE for segment in segments),
    }
    modifiers = [
        modifier
        for enabled, modifier in (
            (result["eci_used"], 1),
            (result["gs1"], 2),
            (result["uri"], 4),
            (result["unicode"], 8),
        )
        if enabled
    ]
    result["symbology_identifier"] = "]h" + (str(modifiers[0]) if modifiers else "0")
    return result


# --------------------------------------------------------------------------
# 矩阵解码: 功能信息 -> 去掩模 -> 栅栏解交织 -> RS 纠错 -> 信息位流
# --------------------------------------------------------------------------
def decode_matrix(modules, dimension=None) -> dict:
    if modules and isinstance(modules[0], (list, tuple)):
        rows = [list(row) for row in modules]
        dimension = len(rows)
        if any(len(row) != dimension for row in rows):
            raise ValueError("矩阵必须是方阵")
        modules = [1 if value else 0 for row in rows for value in row]
    else:
        modules = [1 if value else 0 for value in modules]
        if dimension is None:
            dimension = int(len(modules) ** 0.5)
    if dimension * dimension != len(modules):
        raise ValueError("矩阵维度与模块数不符")

    dimension_version = version_for_dimension(dimension)
    function_information = read_function_information(modules, dimension)
    if function_information["version"] != dimension_version:
        raise FunctionInformationError(
            f"功能信息声明版本 {function_information['version']}, 但矩阵尺寸对应版本 {dimension_version}"
        )

    version = function_information["version"]
    template = build_template(version)
    mask = function_information["mask"]
    unmasked = list(modules)
    for index, flag in enumerate(template.functions):
        if flag:
            continue
        row, column = divmod(index, dimension)
        if mask_bit(mask, row, column):
            unmasked[index] ^= 1

    bit_offsets, _ = data_placement(version)
    total_codewords = len(bit_offsets) // 8
    placed = bytearray(total_codewords)
    for bit_index, offset in enumerate(bit_offsets):
        if unmasked[offset]:
            placed[bit_index >> 3] |= 0x80 >> (bit_index & 7)

    data_codewords, corrected_errors, block_reports = correct_blocks(
        list(placed), version, function_information["level"]
    )
    payload = read_payload(data_codewords)
    payload.update({
        "modules": list(modules),
        "ok": True,
        "format": "han-xin",
        "version": version,
        "dimension": dimension,
        "level": function_information["level"],
        "error_level": EC_LEVELS[function_information["level"]],
        "mask": mask,
        "corrected_errors": corrected_errors,
        "function_information": function_information,
        "rs_blocks": block_reports,
    })
    return payload


# --------------------------------------------------------------------------
# 编码: 文本分段 -> 版本选择 -> RS 纠错 -> 栅栏交织 -> 落模 -> 掩模
# --------------------------------------------------------------------------
def _text_values():
    one, two = {}, {}
    for index in range(10):
        one[str(index)] = index
    for index in range(26):
        one[chr(ord("A") + index)] = 10 + index
    for index in range(26):
        one[chr(ord("a") + index)] = 36 + index
    for value in range(63):
        two[chr(_text_two_byte(value))] = value
    return one, two


TEXT_VALUES_ONE, TEXT_VALUES_TWO = _text_values()


def _character_class(character: str) -> str:
    code = ord(character)
    if 0x30 <= code <= 0x39:
        return "numeric"
    if code < 0x80:
        return "text"
    try:
        encoded = character.encode("gb2312")
    except UnicodeEncodeError:
        encoded = b""
    if len(encoded) == 2 and 0xB0 <= encoded[0] <= 0xF7:
        return "chinese"
    try:
        encoded = character.encode("gb18030")
    except UnicodeEncodeError:
        return "binary"
    if len(encoded) == 2 and 0x81 <= encoded[0] <= 0xFE and 0x40 <= encoded[1] <= 0xFE and encoded[1] != 0x7F:
        return "gb18030"
    return "binary"


def split_runs(text: str):
    runs = []
    for character in text:
        kind = _character_class(character)
        if runs and runs[-1][0] == kind:
            runs[-1][1].append(character)
        else:
            runs.append((kind, [character]))
    return [(kind, "".join(chars)) for kind, chars in runs]


def _write_numeric(writer: BitWriter, digits: str) -> None:
    writer.write(0b0001, 4)
    groups = [int(digits[index:index + 3]) for index in range(0, len(digits), 3)]
    for value in groups:
        writer.write(value, 10)
    final_digits = len(digits) - 3 * (len(groups) - 1)
    writer.write(1020 + final_digits, 10)


def _write_text(writer: BitWriter, chunk: str) -> None:
    writer.write(0b0010, 4)
    submode = 1
    for character in chunk:
        if character in TEXT_VALUES_ONE:
            required, value = 1, TEXT_VALUES_ONE[character]
        elif character in TEXT_VALUES_TWO:
            required, value = 2, TEXT_VALUES_TWO[character]
        else:
            raise HanXinError(f"文本模式无法编码字符 {character!r}")
        if required != submode:
            writer.write(TEXT_SWITCH, 6)
            submode = required
        writer.write(value, 6)
    writer.write(TEXT_TERMINATOR, 6)


def _write_common_chinese(writer: BitWriter, chunk: str) -> None:
    pairs = []
    for character in chunk:
        encoded = character.encode("gb2312")
        first, second = encoded[0], encoded[1]
        if 0xB0 <= first <= 0xD7:
            pairs.append((1, (first - 0xB0) * 94 + (second - 0xA1)))
        else:
            pairs.append((2, (first - 0xD8) * 94 + (second - 0xA1)))
    region = pairs[0][0]
    writer.write(0b0100 if region == 1 else 0b0101, 4)
    for new_region, value in pairs:
        if new_region != region:
            writer.write(CHINESE_SWITCH, 12)
            region = new_region
        writer.write(value, 12)
    writer.write(CHINESE_TERMINATOR, 12)


def _write_gb18030_two_byte(writer: BitWriter, chunk: str) -> None:
    writer.write(0b0110, 4)
    for character in chunk:
        first, second = character.encode("gb18030")
        trail = second - 0x40 if second < 0x7F else second - 0x41
        value = (first - 0x81) * 190 + trail
        if value >= 23940:
            raise HanXinError("GB 18030 双字节码值越界")
        writer.write(value, 15)
    writer.write(GB18030_TWO_BYTE_TERMINATOR, 15)


def _write_binary(writer: BitWriter, data: bytes) -> None:
    if not data:
        return
    writer.write(0b0011, 4)
    writer.write(len(data), 13)
    for value in data:
        writer.write(value, 8)


def build_bitstream(data) -> tuple:
    writer = BitWriter()
    if isinstance(data, (bytes, bytearray, memoryview)):
        data = bytes(data)
        _write_binary(writer, data)
        description = [{"mode": MODE_BINARY, "character_count": len(data)}]
    else:
        description = []
        for kind, chunk in split_runs(data):
            if kind == "numeric":
                _write_numeric(writer, chunk)
            elif kind == "text":
                _write_text(writer, chunk)
            elif kind == "chinese":
                _write_common_chinese(writer, chunk)
            elif kind == "gb18030":
                _write_gb18030_two_byte(writer, chunk)
            else:
                raw = chunk.encode("utf-8")
                _write_binary(writer, raw)
                description.append({"mode": MODE_BINARY, "character_count": len(raw)})
                continue
            description.append({"mode": kind, "character_count": len(chunk)})
    return writer, description


def _mask_penalty(modules, dimension: int) -> int:
    penalty = 0
    for row in range(dimension):
        run_length, previous = 1, modules[row * dimension]
        for column in range(1, dimension):
            value = modules[row * dimension + column]
            if value == previous:
                run_length += 1
            else:
                if run_length >= 5:
                    penalty += 3 + run_length - 5
                run_length, previous = 1, value
        if run_length >= 5:
            penalty += 3 + run_length - 5
    for column in range(dimension):
        run_length, previous = 1, modules[column]
        for row in range(1, dimension):
            value = modules[row * dimension + column]
            if value == previous:
                run_length += 1
            else:
                if run_length >= 5:
                    penalty += 3 + run_length - 5
                run_length, previous = 1, value
        if run_length >= 5:
            penalty += 3 + run_length - 5
    for row in range(dimension - 1):
        for column in range(dimension - 1):
            value = modules[row * dimension + column]
            if (
                value == modules[row * dimension + column + 1]
                and value == modules[(row + 1) * dimension + column]
                and value == modules[(row + 1) * dimension + column + 1]
            ):
                penalty += 3
    dark = sum(modules)
    total = dimension * dimension
    penalty += abs(dark * 100 // total - 50) // 5 * 10
    return penalty


def apply_mask(modules, functions, dimension: int, mask: int):
    result = list(modules)
    for index, flag in enumerate(functions):
        if flag:
            continue
        row, column = divmod(index, dimension)
        if mask_bit(mask, row, column):
            result[index] ^= 1
    return result


def encode(data, error_correction: str = "L1", version=None, mask=None) -> dict:
    """把文本或字节编码成汉信码模块矩阵。"""
    if error_correction not in EC_LEVELS:
        raise ValueError(f"纠错等级必须是 {EC_LEVELS} 之一")
    level = EC_LEVELS.index(error_correction)

    if not data:
        raise ValueError("待编码内容不能为空")
    writer, description = build_bitstream(data)
    needed_bits = len(writer)

    if version is None:
        for candidate in range(MIN_VERSION, MAX_VERSION + 1):
            if rs_structure(candidate, level)["data_bits"] >= needed_bits:
                version = candidate
                break
        else:
            raise HanXinError("内容超出汉信码最大容量")
    elif not MIN_VERSION <= version <= MAX_VERSION:
        raise ValueError(f"版本必须是 {MIN_VERSION}..{MAX_VERSION}")

    structure = rs_structure(version, level)
    if structure["data_bits"] < needed_bits:
        raise HanXinError(
            f"版本 {version} 在 {error_correction} 下只有 {structure['data_bits']} 位, "
            f"需要 {needed_bits} 位"
        )

    data_codewords = bytearray(structure["data_codewords"])
    for index, bit in enumerate(writer.bits):
        if bit:
            data_codewords[index >> 3] |= 0x80 >> (index & 7)

    sequential = bytearray()
    offset = 0
    for batch_index, data_count, correction_count in expand_rs_blocks(version, level):
        block = data_codewords[offset:offset + data_count]
        offset += data_count
        sequential.extend(rs_encode(GF256, block, correction_count, 1))

    placed = to_picket_fence(sequential)
    template = build_template(version)
    dimension = template.dimension
    modules = list(template.modules)
    bit_offsets, _ = data_placement(version)
    for bit_index, target in enumerate(bit_offsets):
        modules[target] = (placed[bit_index >> 3] >> (7 - (bit_index & 7))) & 1

    if mask is None:
        scores = []
        for candidate_mask in range(4):
            masked = apply_mask(modules, template.functions, dimension, candidate_mask)
            scores.append((_mask_penalty(masked, dimension), candidate_mask))
        mask = min(scores)[1]
    elif mask not in (0, 1, 2, 3):
        raise ValueError("掩模号必须是 0..3")

    modules = apply_mask(modules, template.functions, dimension, mask)

    information = encode_function_information(version, level, mask)
    bits = []
    for symbol in information:
        for shift in range(3, -1, -1):
            bits.append((symbol >> shift) & 1)
    bits.extend([0] * FUNCTION_INFO_FILLER_BITS)
    primary, secondary = function_information_coordinates(dimension)
    for coordinates in (primary, secondary):
        for index, (column, row) in enumerate(coordinates):
            modules[row * dimension + column] = bits[index]

    return {
        "modules": modules,
        "dimension": dimension,
        "version": version,
        "level": level,
        "error_level": EC_LEVELS[level],
        "mask": mask,
        "segments": description,
        "data_codewords": structure["data_codewords"],
        "used_bits": needed_bits,
        "capacity_bits": structure["data_bits"],
    }


# --------------------------------------------------------------------------
# 图像输出 (纯标准库写 PNG) 与图像输入 (Pillow 采样)
# --------------------------------------------------------------------------
def _png_chunk(tag: bytes, payload: bytes) -> bytes:
    return (
        struct.pack(">I", len(payload))
        + tag
        + payload
        + struct.pack(">I", zlib.crc32(tag + payload) & 0xFFFFFFFF)
    )


def render_png(modules, dimension: int, scale: int = 10, quiet: int = 2) -> bytes:
    size = (dimension + 2 * quiet) * scale
    raw = bytearray()
    white_row = bytes([0]) + b"\xff" * size
    for y in range(size):
        module_row = y // scale - quiet
        if 0 <= module_row < dimension:
            base = module_row * dimension
            row = bytearray([0])
            for x in range(size):
                module_column = x // scale - quiet
                dark = 0 <= module_column < dimension and modules[base + module_column]
                row.append(0x00 if dark else 0xFF)
            raw += row
        else:
            raw += white_row
    header = struct.pack(">IIBBBBB", size, size, 8, 0, 0, 0, 0)
    return b"".join((
        b"\x89PNG\r\n\x1a\n",
        _png_chunk(b"IHDR", header),
        _png_chunk(b"IDAT", zlib.compress(bytes(raw), 9)),
        _png_chunk(b"IEND", b""),
    ))


def _otsu_threshold(histogram) -> int:
    total = sum(histogram)
    sum_all = sum(level * count for level, count in enumerate(histogram))
    background_weight = 0
    background_sum = 0
    best_threshold, best_variance = 127, -1.0
    for level in range(256):
        background_weight += histogram[level]
        if background_weight == 0:
            continue
        foreground_weight = total - background_weight
        if foreground_weight == 0:
            break
        background_sum += level * histogram[level]
        mean_background = background_sum / background_weight
        mean_foreground = (sum_all - background_sum) / foreground_weight
        variance = background_weight * foreground_weight * (mean_background - mean_foreground) ** 2
        if variance > best_variance:
            best_variance = variance
            best_threshold = level
    return best_threshold


def _sample_grid(mask_image, bbox, dimension: int):
    module = mask_image.crop(bbox).resize((dimension, dimension), _resampling_box())
    return [1 if value < 128 else 0 for value in module.tobytes()]


def _resampling_box():
    from PIL import Image

    resampling = getattr(Image, "Resampling", Image)
    return resampling.BOX


def rotate_grid(modules, dimension: int, quarter_turns: int):
    quarter_turns %= 4
    if quarter_turns == 0:
        return list(modules)
    result = list(modules)
    for _ in range(quarter_turns):
        rotated = [0] * (dimension * dimension)
        for row in range(dimension):
            base = row * dimension
            for column in range(dimension):
                rotated[column * dimension + (dimension - 1 - row)] = result[base + column]
        result = rotated
    return result


def finder_mismatches(modules, dimension: int) -> int:
    mismatches = 0
    for origin_column, origin_row, rotation in (
        (0, 0, 0),
        (dimension - FINDER_SIZE, 0, 1),
        (0, dimension - FINDER_SIZE, 1),
        (dimension - FINDER_SIZE, dimension - FINDER_SIZE, 2),
    ):
        for y in range(FINDER_SIZE):
            for x in range(FINDER_SIZE):
                expected = _finder_value(rotation, x, y)
                actual = modules[(origin_row + y) * dimension + origin_column + x]
                if expected != actual:
                    mismatches += 1
    return mismatches


def image_to_matrix(path: str, forced_dimension=None):
    """把图片采样成模块矩阵, 返回 (modules, dimension, rotation, mismatches)。"""
    try:
        from PIL import Image, ImageOps
    except ImportError as error:  # pragma: no cover
        raise HanXinError("解析图片需要 Pillow: pip install pillow") from error

    with Image.open(path) as handle:
        gray = handle.convert("L")
        gray.load()

    threshold = _otsu_threshold(gray.histogram())
    candidates = []
    for inverted in (False, True):
        mask = gray.point(
            lambda value, inverted=inverted, threshold=threshold: (
                0 if ((value <= threshold) != inverted) else 255
            )
        )
        bbox = ImageOps.invert(mask).getbbox()
        if bbox is None:
            continue
        left, top, right, bottom = bbox
        width, height = right - left, bottom - top
        if width != height or width < 23:
            continue
        dimensions = [forced_dimension] if forced_dimension else [
            dimension_for_version(version) for version in range(MIN_VERSION, MAX_VERSION + 1)
        ]
        for dimension in dimensions:
            module_size = width / dimension
            if not 2 <= module_size <= 256:
                continue
            grid = _sample_grid(mask, bbox, dimension)
            for rotation in range(4):
                rotated = rotate_grid(grid, dimension, rotation)
                candidates.append((finder_mismatches(rotated, dimension), dimension, rotation, rotated))

    if not candidates:
        raise HanXinError("未能从图片中定位汉信码: 请确认图片包含完整的、正方形且高对比的符号")
    candidates.sort(key=lambda item: item[0])
    return candidates


def decode_image(path: str, forced_dimension=None) -> dict:
    candidates = image_to_matrix(path, forced_dimension)
    last_error = None
    for mismatches, dimension, rotation, modules in candidates:
        try:
            result = decode_matrix(modules, dimension)
        except HanXinError as error:
            last_error = error
            continue
        result["finder_mismatches"] = mismatches
        result["rotation"] = rotation * 90
        return result
    raise HanXinError(f"图片中的符号无法解码: {last_error}")


def matrix_from_text(text: str):
    rows = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        rows.append([1 if character in "1#" else 0 for character in line if character in "01#."])
    if not rows:
        raise HanXinError("文本矩阵为空")
    return rows, len(rows)


def matrix_to_text(modules, dimension: int) -> str:
    lines = []
    for row in range(dimension):
        lines.append("".join("██" if modules[row * dimension + column] else "  " for column in range(dimension)))
    return "\n".join(lines)


# --------------------------------------------------------------------------
# 命令行
# --------------------------------------------------------------------------
def _print_result(path: str, result: dict) -> None:
    print(f"文件: {path}")
    if result["text"] is not None:
        print(f"内容: {result['text']}")
    else:
        print(f"内容(hex): {result['bytes'].hex()}")
    print(
        f"版本: {result['version']} ({result['dimension']}x{result['dimension']})  "
        f"纠错: {result['error_level']}  掩模: {result['mask']}  "
        f"纠正码字: {result['corrected_errors']}  旋转: {result.get('rotation', 0)}°  "
        f"符号标识符: {result['symbology_identifier']}"
    )
    modes = ", ".join(
        segment["mode"] + (f"({segment['character_count']})" if segment.get("character_count") else "")
        for segment in result["segments"]
    )
    print(f"数据段: {modes}")


def _selftest() -> int:
    samples = ("HANXIN-SELFTEST", "1234567890", "汉信码 Han Xin 2026")
    failures = 0
    for text in samples:
        result = encode(text)
        decoded = decode_matrix(result["modules"], result["dimension"])
        ok = decoded["text"] == text
        failures += 0 if ok else 1
        print(f"[{'PASS' if ok else 'FAIL'}] 矩阵往返 {result['dimension']}x{result['dimension']} "
              f"{result['error_level']} mask{result['mask']}: {text}")
    try:
        import PIL  # noqa: F401
    except ImportError:
        print("[SKIP] 未安装 Pillow, 跳过图片往返")
        return 1 if failures else 0
    import tempfile
    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, "selftest.png")
        result = encode(samples[2])
        with open(path, "wb") as handle:
            handle.write(render_png(result["modules"], result["dimension"], 8, 2))
        decoded = decode_image(path)
        ok = decoded["text"] == samples[2]
        failures += 0 if ok else 1
        print(f"[{'PASS' if ok else 'FAIL'}] 图片往返 PNG: {decoded['text']}")
    print("自检结束:", "全部通过" if not failures else f"{failures} 项失败")
    return 0 if not failures else 1


def _configure_output() -> None:
    """Windows 控制台默认 GBK, 统一成 UTF-8 并容错, 避免输出非 ASCII 时抛错。"""
    for stream in (sys.stdout, sys.stderr):
        reconfigure = getattr(stream, "reconfigure", None)
        if reconfigure is None:
            continue
        try:
            reconfigure(encoding="utf-8", errors="replace")
        except (ValueError, OSError):
            pass


def main(argv=None) -> int:
    _configure_output()
    parser = argparse.ArgumentParser(
        prog="hanxin.py", description="汉信码 (Han Xin Code) 生成与解析"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    encoder = subparsers.add_parser("encode", help="生成汉信码 PNG")
    encoder.add_argument("text", help="要编码的文本")
    encoder.add_argument("-o", "--output", default="hanxin.png", help="输出 PNG 路径")
    encoder.add_argument("--ec", default="L1", choices=EC_LEVELS, help="纠错等级 (默认 L1)")
    encoder.add_argument("--version", type=int, help="强制版本 1..84 (默认自动选择最小可用版本)")
    encoder.add_argument("--mask", type=int, choices=(0, 1, 2, 3), help="强制掩模号 (默认自动选择)")
    encoder.add_argument("--scale", type=int, default=10, help="每模块像素数 (默认 10)")
    encoder.add_argument("--quiet", type=int, default=2, help="静区宽度(模块, 默认 2)")
    encoder.add_argument("--ascii", action="store_true", help="同时打印模块矩阵")

    subparsers.add_parser("selftest", help="自检: 编码-解码往返 + 图片读写")

    decoder = subparsers.add_parser("decode", help="解析汉信码图片或矩阵文本")
    decoder.add_argument("files", nargs="+", help="图片路径 (PNG/JPEG...) 或 0/1 矩阵文本")
    decoder.add_argument("-j", "--json", action="store_true", help="输出 JSON")
    decoder.add_argument("--ascii", action="store_true", help="同时打印模块矩阵")
    decoder.add_argument("--version", type=int, help="强制版本, 跳过尺寸猜测")

    arguments = parser.parse_args(argv)

    if arguments.command == "selftest":
        return _selftest()

    if arguments.command == "encode":
        try:
            result = encode(
                arguments.text,
                error_correction=arguments.ec,
                version=arguments.version,
                mask=arguments.mask,
            )
            payload = render_png(result["modules"], result["dimension"], arguments.scale, arguments.quiet)
            with open(arguments.output, "wb") as handle:
                handle.write(payload)
        except (HanXinError, ValueError, OSError) as error:
            print(f"生成失败: {error}", file=sys.stderr)
            return 1
        print(
            f"已生成 {arguments.output}  版本 {result['version']} "
            f"({result['dimension']}x{result['dimension']})  纠错 {result['error_level']}  "
            f"掩模 {result['mask']}  使用 {result['used_bits']}/{result['capacity_bits']} 位"
        )
        if arguments.ascii:
            print(matrix_to_text(result["modules"], result["dimension"]))
        return 0

    exit_code = 0
    for path in arguments.files:
        try:
            if path.lower().endswith((".txt", ".mtx")):
                with open(path, "r", encoding="utf-8") as handle:
                    rows, dimension = matrix_from_text(handle.read())
                result = decode_matrix(rows, dimension)
                result.setdefault("rotation", 0)
            else:
                forced = dimension_for_version(arguments.version) if arguments.version else None
                result = decode_image(path, forced)
        except (HanXinError, ValueError, OSError) as error:
            print(f"{path}: 解析失败 - {error}", file=sys.stderr)
            exit_code = 1
            continue
        if arguments.json:
            serializable = dict(result)
            serializable["bytes"] = result["bytes"].hex()
            serializable["segments"] = [
                {key: (value.hex() if isinstance(value, bytes) else value) for key, value in segment.items()}
                for segment in result["segments"]
            ]
            serializable.pop("modules", None)
            print(json.dumps(serializable, ensure_ascii=False, indent=2))
        else:
            _print_result(path, result)
        if arguments.ascii:
            print(matrix_to_text(result["modules"], result["dimension"]))
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
