#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Лабораторная работа №1. Простейшие шифры.
Вариант 5: шифрование и дешифрование методом решётки Кардано.

Ключ - текстовый файл с решёткой 8x8 (8 строк по 8 символов '0'/'1',
1 - отверстие, 0 - нет отверстия); такой файл создаёт режим gen
(вариант 4).

Использование:
    python cardano.py gen     key.txt
    python cardano.py encrypt input.txt key.txt output.txt
    python cardano.py decrypt input.txt key.txt output.txt

Пробелы и символы конца строки не шифруются и остаются на своих местах;
шифруются только остальные символы (буквы).
"""
import random
import sys

N = 8  # размер решётки (чётный)


# ---------------------------------------------------------------- решётка
def rotate(cells):
    """Поворот набора клеток (r, c) на 90 градусов по часовой стрелке."""
    return [(c, N - 1 - r) for r, c in cells]


def generate_grille():
    """Случайная решётка: в каждой из 16 орбит поворота ровно одно отверстие."""
    holes = []
    for r in range(N // 2):
        for c in range(N // 2):
            cell = (r, c)
            for _ in range(random.randrange(4)):
                cell = rotate([cell])[0]
            holes.append(cell)
    grid = [['0'] * N for _ in range(N)]
    for r, c in holes:
        grid[r][c] = '1'
    return [''.join(row) for row in grid]


def load_grille(path):
    with open(path, encoding='utf-8') as f:
        rows = [line.strip() for line in f if line.strip()]
    if len(rows) != N or any(len(r) != N or set(r) - {'0', '1'} for r in rows):
        raise ValueError(f'Решётка должна содержать {N} строк по {N} символов 0/1')
    holes = [(r, c) for r in range(N) for c in range(N) if rows[r][c] == '1']
    # Проверка: при 4 поворотах отверстия покрывают все клетки ровно по разу
    seen = set()
    cur = holes
    for _ in range(4):
        seen.update(cur)
        cur = rotate(cur)
    if len(holes) != N * N // 4 or len(seen) != N * N:
        raise ValueError('Некорректная решётка: отверстия при поворотах '
                         'не покрывают все клетки ровно по одному разу')
    return holes


def cell_order(holes):
    """Порядок заполнения 64 клеток: 4 положения решётки,
    в каждом - отверстия слева направо, сверху вниз."""
    order = []
    cur = sorted(holes)
    for _ in range(4):
        order.extend(cur)
        cur = sorted(rotate(cur))
    return order


# -------------------------------------------------------------- шифрование
def split_text(text):
    """Отделяет шифруемые символы от пробелов и концов строк."""
    return [ch for ch in text if not ch.isspace()]


def merge_text(template, symbols):
    """Вставляет символы обратно в шаблон, сохраняя пробелы и переводы строк."""
    it = iter(symbols)
    return ''.join(ch if ch.isspace() else next(it) for ch in template)


def encrypt_symbols(symbols, order):
    size = N * N
    out = []
    for i in range(0, len(symbols), size):
        block = symbols[i:i + size]
        grid = {}
        for ch, pos in zip(block, order):  # неполный блок заполняет первые клетки
            grid[pos] = ch
        for r in range(N):
            for c in range(N):
                if (r, c) in grid:
                    out.append(grid[(r, c)])
    return out


def decrypt_symbols(symbols, order):
    size = N * N
    out = []
    for i in range(0, len(symbols), size):
        block = symbols[i:i + size]
        used = sorted(order[:len(block)])      # занятые клетки в порядке чтения
        grid = dict(zip(used, block))
        out.extend(grid[pos] for pos in order[:len(block)])
    return out


# ------------------------------------------------------------------- main
def main(argv):
    if len(argv) >= 3 and argv[1] == 'gen':
        with open(argv[2], 'w', encoding='utf-8') as f:
            f.write('\n'.join(generate_grille()) + '\n')
        print('Решётка записана в', argv[2])
    elif len(argv) == 5 and argv[1] in ('encrypt', 'decrypt'):
        with open(argv[2], encoding='utf-8') as f:
            text = f.read()
        order = cell_order(load_grille(argv[3]))
        func = encrypt_symbols if argv[1] == 'encrypt' else decrypt_symbols
        result = merge_text(text, func(split_text(text), order))
        with open(argv[4], 'w', encoding='utf-8') as f:
            f.write(result)
        print('Готово:', argv[4])
    else:
        print(__doc__)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
