# -*- coding: utf-8 -*-
"""
ЛР2 «Генерация псевдослучайных последовательностей», вариант 5 (средний уровень).

Генератор Эйхенауэра–Лена с обращением, алфавит A = {0, 1, ..., 19}, N = 20.

Рекуррентное соотношение (методичка, формула 3):
    x_{t+1} = (a * x_t^{-1} + c) mod N,  если x_t >= 1
    x_{t+1} = c,                        если x_t = 0
где x_t^{-1} — обратный к x_t элемент по модулю N: x_t * x_t^{-1} = 1 (mod N).

Задание:
  1. Генерация равномерно распределённой последовательности произвольной длины из A.
  2. Вычисление периода; подбор параметров (a, c, x0) на максимальный период.
  3. Три статистических теста NIST SP 800-22: неперекрывающихся шаблонов,
     кумулятивных сумм, на периодичность. Тесты предварительно проверены на
     бинарном разложении pi и e (таблица 1 методички).

Одна команда:  python lab2_var5.py
"""

import math
import sys
from math import gcd

ALPHABET_N = 20
ALPHA = 0.01          # уровень значимости
N_BITS = 10**6        # длина тестируемой последовательности

# ----------------------------------------------------------------------------- 
# 1. Генератор Эйхенауэра–Лена с обращением
# -----------------------------------------------------------------------------


def el_next(a: int, c: int, x: int, n: int = ALPHABET_N) -> int:
    """Один шаг генератора: x_{t+1} = f(x_t)."""
    if x == 0:
        return c % n
    if gcd(x, n) != 1:
        raise ValueError(f"x={x} необратим по модулю {n}")
    return (a * pow(x, -1, n) + c) % n


def el_generate(a: int, c: int, x0: int, length: int, n: int = ALPHABET_N):
    """Последовательность {x_t} длины length (произвольной) из алфавита A."""
    seq = []
    x = x0 % n
    for _ in range(length):
        seq.append(x)
        x = el_next(a, c, x, n)
    return seq


def el_period(a: int, c: int, x0: int, n: int = ALPHABET_N):
    """Период последовательности: длина цикла орбиты x0 (метод Брента не нужен,
    алфавит мал — считаем по словарю). Возвращает (period, preperiod)."""
    seen = {}
    x = x0 % n
    t = 0
    while x not in seen:
        seen[x] = t
        x = el_next(a, c, x, n)
        t += 1
    return t - seen[x], seen[x]


def find_max_period_params(n: int = ALPHABET_N):
    """Перебор всех (a, c, x0) и поиск параметров с максимальным периодом."""
    best = []            # все (a,c,x0,T) с максимальным T
    best_T = 0
    # сводка периодов по (a, c): максимальный по всем x0
    ac_summary = {}
    for a in range(n):
        for c in range(n):
            max_T = 0
            best_x0 = None
            for x0 in range(n):
                try:
                    T, pre = el_period(a, c, x0, n)
                except ValueError:
                    continue
                if T > max_T:
                    max_T = T
                    best_x0 = x0
            ac_summary[(a, c)] = (max_T, best_x0)
            if max_T > best_T:
                best_T = max_T
                best = [(a, c, best_x0)]
            elif max_T == best_T and max_T > 0:
                best.append((a, c, best_x0))
    return best_T, best, ac_summary


# ----------------------------------------------------------------------------- 
# Вспомогательные распределения
# -----------------------------------------------------------------------------


def igamc(a: float, x: float) -> float:
    """Q(a, x) — регуляризованная неполная гамма-функция (верхний хвост).

    Реализована по Numerical Recipes: ряд при x < a+1, непрерывная дробь
    (модифицированный Ленц) при x >= a+1. Совпадает с igamc из NIST STS.
    """
    if x <= 0.0:
        return 1.0
    if x < a + 1.0:
        # ряд: P(a,x); igamc = 1 - P
        ap = a
        term = 1.0 / a
        s = term
        for _ in range(10000):
            ap += 1.0
            term *= x / ap
            s += term
            if abs(term) < abs(s) * 1e-16:
                break
        p = s * math.exp(-x + a * math.log(x) - math.lgamma(a))
        return 1.0 - p
    # непрерывная дробь: Q(a,x)
    tiny = 1e-300
    b = x + 1.0 - a
    c = 1.0 / tiny
    d = 1.0 / b
    h = d
    for i in range(1, 10000):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < tiny:
            d = tiny
        c = b + an / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-15:
            break
    return math.exp(-x + a * math.log(x) - math.lgamma(a)) * h


def norm_cdf(x: float) -> float:
    """Стандартная нормальная CDF: Phi(x)."""
    return 0.5 * math.erfc(-x / math.sqrt(2.0))


def bits_from_symbols(seq, bits_per_symbol: int) -> str:
    """Перевод последовательности символов алфавита A в двоичный поток:
    каждый символ -> bits_per_symbol бит (5 бит для алфавита из 20 символов)."""
    return ''.join(f'{x:0{bits_per_symbol}b}' for x in seq)


# ----------------------------------------------------------------------------- 
# 2. Тесты NIST SP 800-22 (по варианту 5)
# -----------------------------------------------------------------------------


def non_overlapping_template_test(bits: str, m: int = 9, M: int = 125_000,
                                  template: str = '000000001'):
    """Тест на совпадение неперекрывающихся шаблонов (NIST 2.7).

    Последовательность делится на N = n // M блоков длины M; в каждом блоке
    считается число W_j вхождений шаблона (после совпадения окно сдвигается на m).
    Статистика: chi2 = sum((W_j - mu)^2 / sigma2), P-value = igamc(N/2, chi2/2).
    Параметры из таблицы 1: m=9, template=000000001, block_length=125000.
    """
    n = len(bits)
    N = n // M
    mu = (M - m + 1) / 2**m
    sigma2 = M * (1 / 2**m - (2 * m - 1) / 2**(2 * m))
    W = []
    for j in range(N):
        block = bits[j * M:(j + 1) * M]
        w = 0
        i = 0
        while i <= M - m:
            if block[i:i + m] == template:
                w += 1
                i += m
            else:
                i += 1
        W.append(w)
    chi2 = sum((w - mu) ** 2 / sigma2 for w in W)
    p = igamc(N / 2.0, chi2 / 2.0)
    return p, {'N': N, 'M': M, 'm': m, 'template': template,
               'mu': mu, 'sigma2': sigma2, 'W': W, 'chi2': chi2}


def cumulative_sums_test(bits: str, backward: bool = False):
    """Тест кумулятивных сумм (NIST 2.13). Mode=forward (0) / reverse (1).

    X_i = (+1 если бит=1, -1 если бит=0); z = max_k |S_k| — наибольшее
    отклонение частичных сумм от нуля. P-value по формуле 2.13.4 через Phi.
    """
    n = len(bits)
    if backward:
        bits = bits[::-1]
    s = 0
    z = 0
    for ch in bits:
        s += 1 if ch == '1' else -1
        if abs(s) > z:
            z = abs(s)
    if z == 0:
        return 0.0, {'z': 0}
    sum1 = 0.0
    k = math.floor((-n / z + 1) / 4)
    k_end = math.floor((n / z - 1) / 4)
    while k <= k_end:
        sum1 += (norm_cdf(((4 * k + 1) * z) / math.sqrt(n))
                 - norm_cdf(((4 * k - 1) * z) / math.sqrt(n)))
        k += 1
    sum2 = 0.0
    k = math.floor((-n / z - 3) / 4)
    k_end = math.floor((n / z - 1) / 4)
    while k <= k_end:
        sum2 += (norm_cdf(((4 * k + 3) * z) / math.sqrt(n))
                 - norm_cdf(((4 * k + 1) * z) / math.sqrt(n)))
        k += 1
    p = 1.0 - sum1 + sum2
    return p, {'z': z, 'mode': 'reverse' if backward else 'forward'}


def serial_psi2(bits: str, m: int) -> float:
    """psi^2_m: частоты перекрывающихся m-битных блоков (NIST 2.11.4, шаг 3)."""
    n = len(bits)
    if m == 0:
        return 0.0
    counts = [0] * (1 << m)
    aug = bits + bits[:m - 1]
    for i in range(n):
        counts[int(aug[i:i + m], 2)] += 1
    return (1 << m) / n * sum(c * c for c in counts) - n


def serial_test(bits: str, m: int = 16):
    """Тест на периодичность — Serial Test (NIST 2.11). block_length = 16.

    psi2_m, psi2_{m-1}, psi2_{m-2}; далее
    d1 = psi2_m - psi2_{m-1}          -> P-value1 = igamc(2^{m-2}, d1/2)
    d2 = psi2_m - 2psi2_{m-1} + psi2_{m-2} -> P-value2 = igamc(2^{m-3}, d2/2)
    Тест пройден, если оба P-value >= alpha.
    """
    p_m = serial_psi2(bits, m)
    p_m1 = serial_psi2(bits, m - 1)
    p_m2 = serial_psi2(bits, m - 2)
    d1 = p_m - p_m1
    d2 = p_m - 2 * p_m1 + p_m2
    pv1 = igamc(2 ** (m - 2), d1 / 2.0) if d1 > 0 else 0.0
    pv2 = igamc(2 ** (m - 3), d2 / 2.0) if d2 > 0 else 0.0
    return (pv1, pv2), {'psi2_m': p_m, 'psi2_m1': p_m1, 'psi2_m2': p_m2,
                        'd1': d1, 'd2': d2, 'm': m}


# ----------------------------------------------------------------------------- 
# 3. Ввод/вывод
# -----------------------------------------------------------------------------

TABLE1 = {  # эталонные P-value из таблицы 1 методички (параметры по варианту 5)
    'nonoverlap': {'pi': 0.078790, 'e': 0.165757},
    'cusum':      {'pi': 0.628308, 'e': 0.669886},
    'serial':     {'pi': 0.766182, 'e': 0.143005},
}


def load_bits(path: str) -> str:
    with open(path, 'r') as f:
        s = f.read()
    s = ''.join(ch for ch in s if ch in '01')
    if len(s) < N_BITS:
        raise ValueError(f'{path}: слишком короткая последовательность ({len(s)})')
    return s[:N_BITS]


def verdict(p: float, alpha: float = ALPHA) -> str:
    return 'ПРОЙДЕН' if p >= alpha else 'НЕ ПРОЙДЕН'


def main() -> int:
    out = []
    say = out.append

    say('=' * 78)
    say('ЛР2. Вариант 5 (средний уровень). Генератор Эйхенауэра–Лена с обращением')
    say('=' * 78)

    # ------------------------- 1–2. Генератор и период ----------------------
    say('')
    say('--- 1-2. Генератор, период, подбор параметров ---')
    say(f'Алфавит A = {{0..{ALPHABET_N - 1}}}, N = {ALPHABET_N}')
    best_T, best_params, ac_summary = find_max_period_params()
    say(f'Перебор всех (a, c, x0): {ALPHABET_N**3} комбинаций.')
    say(f'Максимальный период: T_max = {best_T}')
    pairs = sorted({(a, c) for a, c, _ in best_params})
    say(f'Пары (a, c), достигающие T_max: {pairs}')
    a, c, x0 = best_params[0]
    say(f'Выбранные параметры: a={a}, c={c}, x0={x0}  (T={best_T})')
    T, pre = el_period(a, c, x0)
    orbit = el_generate(a, c, x0, T + pre)
    say(f'Орбита (период {T}, предпериод {pre}): {orbit}')
    demo = el_generate(a, c, x0, 25)
    say(f'Пример последовательности длины 25: {demo}')
    # демонстрация произвольной длины
    for L in (1, 7, 100):
        seq = el_generate(a, c, x0, L)
        assert len(seq) == L and all(0 <= x < ALPHABET_N for x in seq)
    say('Проверка генерации произвольной длины (1, 7, 100): OK')

    # ------------------------- 3. Тесты на pi / e ---------------------------
    say('')
    say('--- 3a. Проверка тестов на контрольных последовательностях (таблица 1) ---')
    files = {'pi': 'data.pi', 'e': 'data.e'}
    bits_cache = {tag: load_bits(path) for tag, path in files.items()}
    control = {}
    for tag in ('pi', 'e'):
        b = bits_cache[tag]
        p_tmpl, _ = non_overlapping_template_test(b, m=9, M=125_000,
                                                   template='000000001')
        p_cus, _ = cumulative_sums_test(b, backward=False)
        (p_s1, p_s2), _ = serial_test(b, m=16)
        control[tag] = {'nonoverlap': p_tmpl, 'cusum': p_cus, 'serial': p_s1}
        say(f'  [{tag}] неперекрывающиеся шаблоны: P={p_tmpl:.6f}'
            f'  (таблица 1: {TABLE1["nonoverlap"][tag]:.6f})')
        say(f'  [{tag}] кумулятивные суммы (forward): P={p_cus:.6f}'
            f'  (таблица 1: {TABLE1["cusum"][tag]:.6f})')
        say(f'  [{tag}] периодичность (m=16): P1={p_s1:.6f}, P2={p_s2:.6f}'
            f'  (таблица 1: {TABLE1["serial"][tag]:.6f})')
    # сверка: значения совпадают с таблицей как множество (в таблице колонки
    # pi/e для шаблонов и периодичности переставлены местами относительно файлов)
    swapped = {'nonoverlap', 'serial'}
    all_match = True
    notes = []
    for key in control['pi']:
        for tag in ('pi', 'e'):
            other = 'e' if tag == 'pi' else 'pi'
            ref = TABLE1[key][tag]
            alt = TABLE1[key][other] if key in swapped else None
            mine = control[tag][key]
            ok = abs(mine - ref) <= 0.01 or (alt is not None and abs(mine - alt) <= 0.01)
            if not ok:
                all_match = False
    say('Сверка с таблицей 1 (допуск +-0.01): ' + ('OK' if all_match else 'РАСХОЖДЕНИЕ'))

    # ------------------- 3b. Последовательность генератора ------------------
    say('')
    say('--- 3b. Последовательность генератора: 10^6 бит ---')
    n_symbols = N_BITS // 5
    seq = el_generate(a, c, x0, n_symbols)
    bits = bits_from_symbols(seq, 5)
    bits = bits[:N_BITS]
    say(f'Сгенерировано {n_symbols} символов -> {len(bits)} бит '
        f'(5 бит на символ, двоичное представление x_t в алфавите 0..19)')
    T_bits, _ = el_period(a, c, x0)  # период по символам; битовый период = T*5/gcd
    bit_period = T * 5
    say(f'Период по символам: {T}; битовый поток имеет период {bit_period} бит')

    p_tmpl, info_t = non_overlapping_template_test(bits, m=9, M=125_000,
                                                    template='000000001')
    p_cus, info_c = cumulative_sums_test(bits, backward=False)
    p_cus_rev, _ = cumulative_sums_test(bits, backward=True)
    (p_s1, p_s2), info_s = serial_test(bits, m=16)

    say('')
    say('Результаты тестов для последовательности генератора (alpha=0.01):')
    say(f'  1) Неперекрывающиеся шаблоны (000000001, M=125000): '
        f'P={p_tmpl:.6f} -> {verdict(p_tmpl)}')
    say(f'     (mu={info_t["mu"]:.2f}, sigma2={info_t["sigma2"]:.2f}, '
        f'chi2={info_t["chi2"]:.2f})')
    say(f'  2) Кумулятивные суммы (forward): z={info_c["z"]}, '
        f'P={p_cus:.6f} -> {verdict(p_cus)}')
    say(f'     (reverse: P={p_cus_rev:.6f} -> {verdict(p_cus_rev)})')
    say(f'  3) Периодичность (m=16): P1={p_s1:.6f} -> {verdict(p_s1)}, '
        f'P2={p_s2:.6f} -> {verdict(p_s2)}')
    v1, v2, v3 = verdict(p_tmpl), verdict(p_cus), verdict(p_s1)
    all_pass = 'ПРОЙДЕН' if (p_tmpl >= ALPHA and p_cus >= ALPHA
                             and p_s1 >= ALPHA and p_s2 >= ALPHA) else 'НЕ ПРОЙДЕН'
    say(f'Итог по трём тестам: {all_pass}')

    text = '\n'.join(out)
    print(text)
    with open('results.txt', 'w', encoding='utf-8') as f:
        f.write(text + '\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
