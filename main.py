import math
import random
from sympy import isprime
import numpy as np


def inverse(a, n):
    return pow(a, -1, n)


class Point:
    """
    Point:
        Координаты точки p=(x,y).
    """

    def __init__(self, x: int, y: int) -> None:
        self.x = x
        self.y = y

    def __eq__(self, other) -> bool:
        return type(self) == type(other) and (self.x == other.x) and (self.y == other.y)

    def __repr__(self):
        return f"(x,y) = ({self.x}, {self.y})"


class EllipticCurve:
    '''
    Эллиптическая кривая вида y^2 = x^3 + ax + b над Z/NZ
    '''

    def __init__(self, a: int, b: int, N: int) -> None:
        self.a = a
        self.b = b
        self.N = N
        self.div = 1

    def sum(self, p: Point, q: Point) -> Point:
        if type(None) == type(p):
            return q
        if type(None) == type(q):
            return p

        if p == q:
            # При делении на 0 в Z/NZ возвращает None
            if math.gcd(2 * p.y, self.N) != 1:
                self.div = math.gcd(2 * p.y, self.N)
                return None

            slope = (3 * p.x * p.x + self.a) * inverse(2 * p.y, self.N)
            x = (slope * slope - 2 * p.x) % self.N
            y = (slope * (p.x - x) - p.y) % self.N
            return Point(x, y)

        else:
            # При делении на 0 в Z/NZ возвращает None
            if math.gcd(q.x - p.x, self.N) != 1:
                self.div = math.gcd(q.x - p.x, self.N)
                return None

            slope = (q.y - p.y) * inverse(q.x - p.x, self.N)
            x3 = (slope * slope - p.x - q.x) % self.N
            y3 = (slope * (p.x - x3) - p.y) % self.N
            return Point(x3, y3)

    def mult(self, p: Point, k: int) -> Point:
        if p == None:
            return None

        binary_expansion = bin(k)[2:]
        m = len(binary_expansion)

        powers_of_two = [p]
        for _ in range(m - 1):
            new_point = self.sum(powers_of_two[-1], powers_of_two[-1])
            powers_of_two.append(new_point)

        q = None
        for i, v in enumerate(binary_expansion):
            if v == "1":
                q = self.sum(powers_of_two[i], q)
                if q == None:
                    return None
        return q


def lenstra(N):
    """
    Алгоритм Ленстры для факторизации. Написан по конспекту
    https://drive.google.com/drive/folders/1TkqHFzS6HrKXgrXgKE4BVkcSZurZE3Q6 [2022.12.10]
    :param N: Факторизируемое число
    :return: Нетривиальный делитель
    """
    # 1. [Порождение кривой]

    # выбираем случайные A, a, b
    a, b, A = [random.randint(1, N) for _ in range(3)]
    # print(f'\na {a}, b {b}, A {A}')

    # P = (a, b)
    point = Point(a, b)
    # print(f'\npoint.x {point.x}, point.y {point.y}')

    # B = b^2 - a^3 - A*a (mod N)
    B = (b * b - a * a * a - A * a) % N
    # print(f'\nb {B}')

    # E: y^2 = x^3 + Ax + B - порожденная кривая
    E = EllipticCurve(A, B, N)

    # print(f'\nE.a {E.a}, E.b {E.b}, E.N {E.N}, E.div {E.div}')

    # 2. [Цикл]

    '''
    for 2 <= j <= L
        P = jP

    if jP == None:
        E.div = найденный элемент
        if E.div == N
            следующая итерация
        else
        E.div | N - нетривиальный делитель
        return E.div
    '''

    i = 2
    while E.div == 1:
        if point == None:
            break
        point = E.mult(point, i)
        # print(f'E.div {E.div}')
        # print(f'\npoint.x {point.x}, point.y {point.y}')
        i += 1
    return E.div


def get_next(N):
    """
    Вспомогательная функция, проверяет, нужно ли факторизовать число дальше
    """
    # число является простым, просто возвращаем
    if isprime(N):
        return [N]
    # число является составным, нужно факторизовать
    # print(f'\nN - {N}')
    return factorization(N)


def factorization(N):
    """
    Факторизирует число алгоритмом Ленстры
    :param N: Факторизируемое число
    :return: лист простых множителей
    """
    factorizations = []

    if isprime(N):
        factorizations.append(N)
        return N

    # получили нетривиальный делитель, возможно простой
    div = lenstra(N)

    # если нетривиальный делитель равен числу
    while div == N:
        div = lenstra(N)

    # print(f'\nlenstra {div}')

    # print(f'lenstra {div}')
    # нужно факторизовать полученное число
    factorizations += get_next(div)

    # так же есть число N / lenstra, которое также является его делителем
    # если оно простое, get_next вернет его же. Если оно равно 1, то лучше скипнуть

    numb = N // div
    # print(f'\nnumb {numb}, N {N}, numbb {numbb}')

    if not numb == 1:
        factorizations += get_next(numb)

    return factorizations


def to_dict(data):
    """
    Конвертирует лист в словарь по принципу: ключ, количество повторений
    """

    dict_data = {key: data.count(key) for key in set(data)}
    return dict(sorted(dict_data.items()))


if __name__ == "__main__":
    N = 1234569454241865415611952418521
    if not isprime(N):
        print(f'{to_dict(factorization(N))}')
