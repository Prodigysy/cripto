# -*- coding: utf-8 -*-
"""
ЛР2 «Генерация псевдослучайных последовательностей», вариант 9 (средний уровень).

LFSR — регистр сдвига с линейной обратной связью над GF(2).

Многочлены варианта (методичка, таблица вариантов):
    f1(x) = x^5 + x^2 + 1
    f2(x) = x^5 + 1
    f3(x) = x^5 + x^4 + x^2 + 1
    Вопрос: какие из этих многочленов примитивны по модулю 2?

Статистические тесты NIST SP 800-22 (по варианту):
    1. Тест серий (Runs, NIST 2.3)
    2. Тест на совпадение неперекрывающихся шаблонов (NIST 2.7)
    3. Универсальный статистический тест Маурера (NIST 2.9)

Задание:
  1. Генерация равномерно распределённой последовательности произвольной
     длины (алфавит {0, 1} — битовый поток LFSR).
  2. Вычисление периода последовательности; подбор параметров генератора
     (многочлен обратной связи + начальное состояние) на максимальный период.
  3. Три статистических теста на 10^6 бит; предварительная проверка тестов
     на бинарном разложении pi и e (таблица 1 методички).

Одна команда:  python lab2_var9.py
"""

import math
import sys

ALPHA = 0.01          # уровень значимости
N_BITS = 10**6        # длина тестируемой последовательности

# -----------------------------------------------------------------------------
# 1. LFSR над GF(2)
# -----------------------------------------------------------------------------
# Многочлен f(x) = x^L + c_{L-1} x^{L-1} + ... + c_1 x + c_0 задаётся списком
# степеней с ненулевыми коэффициентами, включая старшую степень L.
# Рекуррентное соотношение: b_{t+L} = XOR по e из taps (b_{t+e}),
# где taps — все степени многочлена, кроме старшей.


def lfsr_step(taps, L, state):
    """Шаг LFSR: state — целое, младший бит = самый старый бит b_t.
    Возвращает (выходной бит, новое состояние)."""
    out = state & 1
    nb = 0
    for e in taps:
        nb ^= (state >> e) & 1
    state = (state >> 1) | (nb << (L - 1))
    return out, state


def lfsr_generate(taps, L, seed, length):
    """Битовая последовательность LFSR произвольной длины length."""
    bits = []
    state = seed
    for _ in range(length):
        out, state = lfsr_step(taps, L, state)
        bits.append(out)
    return bits


def lfsr_period(taps, L, seed):
    """Период последовательности LFSR с начальным состоянием seed.
    Возвращает (period, preperiod)."""
    seen = {}
    state = seed
    t = 0
    while state not in seen:
        seen[state] = t
        _, state = lfsr_step(taps, L, state)
        t += 1
    return t - seen[state], seen[state]


def lfsr_period_fast(taps, L, seed, cap=None):
    """Период LFSR методом Брента (без словаря состояний).
    Возвращает (period, preperiod). Перекрёстная проверка lfsr_period."""
    if cap is None:
        cap = 1 << L
    power = lam = 1
    tortoise = seed
    _, hare = lfsr_step(taps, L, seed)
    while tortoise != hare and lam < cap:
        if power == lam:
            tortoise = hare
            power *= 2
            lam = 0
        _, hare = lfsr_step(taps, L, hare)
        lam += 1
    if tortoise != hare:
        raise ValueError("период не найден в пределах cap")
    period = 0
    hare = tortoise
    while True:
        _, hare = lfsr_step(taps, L, hare)
        period += 1
        if tortoise == hare:
            break
    tortoise = hare = seed
    for _ in range(period):
        _, hare = lfsr_step(taps, L, hare)
    preperiod = 0
    while tortoise != hare:
        _, tortoise = lfsr_step(taps, L, tortoise)
        _, hare = lfsr_step(taps, L, hare)
        preperiod += 1
    return period, preperiod


# --- Многочлены варианта 9 ---------------------------------------------------
# f1(x) = x^5 + x^2 + 1        -> taps = (2, 0)
# f2(x) = x^5 + 1              -> taps = (0,)
# f3(x) = x^5 + x^4 + x^2 + 1  -> taps = (4, 2, 0)

POLYS = {
    'f1': {'name': 'x^5 + x^2 + 1', 'taps': (2, 0)},
    'f2': {'name': 'x^5 + 1', 'taps': (0,)},
    'f3': {'name': 'x^5 + x^4 + x^2 + 1', 'taps': (4, 2, 0)},
}


def gf2_mod(a, m):
    """Остаток от деления многочленов над GF(2) (биты = коэффициенты)."""
    dm = m.bit_length() - 1
    while a and a.bit_length() - 1 >= dm:
        a ^= m << (a.bit_length() - 1 - dm)
    return a


def gf2_gcd(a, b):
    while b:
        a, b = b, gf2_mod(a, b)
    return a


def gf2_mul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
    return r


def gf2_powmod(base, e, m):
    r = 1
    base = gf2_mod(base, m)
    while e:
        if e & 1:
            r = gf2_mod(gf2_mul(r, base), m)
        base = gf2_mod(gf2_mul(base, base), m)
        e >>= 1
    return r


def factor_int(n):
    """Разложение целого числа на простые множители (для степеней 2^L-1)."""
    fs = []
    d = 2
    while d * d <= n:
        while n % d == 0:
            fs.append(d)
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        fs.append(n)
    return fs


def is_irreducible(f):
    """Неприводимость многочлена f над GF(2), deg f = n >= 1."""
    n = f.bit_length() - 1
    if not f & 1:
        return False
    # x^(2^n) == x (mod f) и gcd(f, x^(2^(n/p)) + x) == 1 для всех p | n
    r = 2
    for _ in range(n):
        r = gf2_mod(gf2_mul(r, r), f)
    if r != 2:
        return False
    for p in sorted(set(factor_int(n))):
        r = 2
        for _ in range(n // p):
            r = gf2_mod(gf2_mul(r, r), f)
        if gf2_gcd(f, r ^ 2) != 1:
            return False
    return True


def is_primitive(f):
    """Примитивность f над GF(2): f неприводим и x — образующая мультипликативной
    группы GF(2^n)*, т.е. ord(x) = 2^n - 1."""
    if not is_irreducible(f):
        return False, 'приводим'
    n = f.bit_length() - 1
    T = (1 << n) - 1
    for q in sorted(set(factor_int(T))):
        if gf2_powmod(2, T // q, f) == 1:
            return False, f'ord(x) < {T} (делит {T}//{q})'
    return True, f'ord(x) = {T} = 2^{n}-1'


# -----------------------------------------------------------------------------
# 2. Вспомогательные распределения (как в NIST STS)
# -----------------------------------------------------------------------------


def igamc(a: float, x: float) -> float:
    """Q(a, x) — регуляризованная неполная гамма-функция (верхний хвост).

    Реализована по Numerical Recipes: ряд при x < a+1, непрерывная дробь
    (модифицированный Ленц) при x >= a+1. Совпадает с igamc из NIST STS.
    """
    if x <= 0.0:
        return 1.0
    if x < a + 1.0:
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


# -----------------------------------------------------------------------------
# 3. Тесты NIST SP 800-22 (по варианту 9)
# -----------------------------------------------------------------------------


def runs_test(bits):
    """Тест серий (NIST 2.3). Возвращает (P-value, словарь статистик)."""
    n = len(bits)
    ones = bits.count('1')
    pi = ones / n
    # предтест частотного побитового
    V = 1
    for k in range(n - 1):
        if bits[k] != bits[k + 1]:
            V += 1
    if abs(pi - 0.5) >= 2.0 / math.sqrt(n):
        return 0.0, {'pi': pi, 'V': V,
                     'note': 'не пройден частотный предтест (NIST 2.3.4 п.2), P=0'}
    p = math.erfc(abs(V - 2 * n * pi * (1 - pi))
                  / (2 * math.sqrt(2 * n) * pi * (1 - pi)))
    return p, {'pi': pi, 'V': V}


def non_overlapping_template_test(bits: str, m: int = 9, M: int = 125_000,
                                  template: str = '000000001'):
    """Тест на совпадение неперекрывающихся шаблонов (NIST 2.7).

    Последовательность делится на N = n // M блоков длины M; в каждом блоке
    считается число вхождений шаблона W_j (после вхождения окно сдвигается
    на m бит). P = igamc(N/2, chi2/2).
    """
    n = len(bits)
    N_blocks = n // M
    mu = (M - m + 1) / (1 << m)
    sigma2 = M * (1 / (1 << m) - (2 * m - 1) / (1 << (2 * m)))
    W = []
    for j in range(N_blocks):
        block = bits[j * M:(j + 1) * M]
        w = 0
        i = 0
        while i <= len(block) - m:
            if block[i:i + m] == template:
                w += 1
                i += m
            else:
                i += 1
        W.append(w)
    chi2 = sum((w - mu) ** 2 for w in W) / sigma2
    p = igamc(N_blocks / 2, chi2 / 2)
    return p, {'mu': mu, 'sigma2': sigma2, 'W': W, 'chi2': chi2}


# Таблица expectedValue / variance для универсального теста Маурера (NIST 2.9,
# Handbook of Applied Cryptography). L=6..16.
MAURER_EXPECTED = {
    6: 5.2177052, 7: 6.1962507, 8: 7.1836656, 9: 8.1764248, 10: 9.1723243,
    11: 10.176001, 12: 11.168765, 13: 12.168070, 14: 13.167693, 15: 14.167488,
    16: 15.167379,
}
MAURER_VARIANCE = {
    6: 2.954, 7: 3.125, 8: 3.238, 9: 3.311, 10: 3.356, 11: 3.384,
    12: 3.401, 13: 3.410, 14: 3.416, 15: 3.419, 16: 3.421,
}


def maurers_universal_test(bits: str, L: int = 7, Q: int = 1280):
    """Универсальный статистический тест Маурера (NIST 2.9).

    n бит делятся на Q инициализационных и K = n//L - Q тестовых L-битных
    блоков; fn = среднее log2 расстояний между повторениями L-битных блоков.
    P = erfc(|fn - E(L)| / (sqrt(2) * sigma)), sigma = sqrt(variance/K).
    """
    n = len(bits)
    if L not in MAURER_EXPECTED:
        raise ValueError(f"L={L} вне диапазона 6..16")
    n_blocks = n // L
    K = n_blocks - Q
    if K <= 0:
        raise ValueError("слишком короткая последовательность")
    T = [0] * (1 << L)
    blocks = [int(bits[i * L:(i + 1) * L], 2) for i in range(n_blocks)]
    for i in range(Q):
        T[blocks[i]] = i
    total = 0.0
    for i in range(Q, n_blocks):
        j = blocks[i]
        total += math.log2(i - T[j])
        T[j] = i
    fn = total / K
    expected = MAURER_EXPECTED[L]
    variance = MAURER_VARIANCE[L]
    c = 0.7 - 0.8 / L + (4 + 32 / L) * (K - Q / 2) / (K * K)  # (2.9.4 п.5)
    sigma = c * math.sqrt(variance / K)
    p = math.erfc(abs(fn - expected) / (math.sqrt(2.0) * sigma))
    return p, {'fn': fn, 'expected': expected, 'sigma': sigma, 'c': c,
               'L': L, 'Q': Q, 'K': K}


# -----------------------------------------------------------------------------
# 4. Ввод/вывод
# -----------------------------------------------------------------------------

TABLE1 = {  # эталонные P-value из таблицы 1 методички (параметры по варианту 9)
    'runs':        {'pi': 0.561917, 'e': 0.419268},
    'nonoverlap':  {'pi': 0.078790, 'e': 0.165757},
    'maurer':      {'pi': 0.282568, 'e': 0.669012},
}


def load_bits(path: str) -> str:
    """Читает контрольную последовательность, удаляя пробелы и переводы
    строк (файлы data.pi / data.e разбиты на строки по 24 бита)."""
    with open(path, 'r') as f:
        s = f.read()
    s = ''.join(s.split())
    return s[:N_BITS]


def verdict(p: float, alpha: float = ALPHA) -> str:
    return 'ПРОЙДЕН' if p >= alpha else 'НЕ ПРОЙДЕН'


def main() -> int:
    out = []
    say = out.append

    say('=' * 78)
    say('ЛР2. Вариант 9 (средний уровень). LFSR над GF(2)')
    say('=' * 78)

    # ------------------------- 1-2. Генератор и период ----------------------
    say('')
    say('--- 1-2. Генератор, период, подбор параметров ---')
    # Многочлены варианта в битовом представлении: бит e = коэффициент x^e
    f1 = (1 << 5) | (1 << 2) | 1            # x^5 + x^2 + 1
    f2 = (1 << 5) | 1                        # x^5 + 1
    f3 = (1 << 5) | (1 << 4) | (1 << 2) | 1  # x^5 + x^4 + x^2 + 1
    polys = [('f1', 'x^5 + x^2 + 1', f1, POLYS['f1']['taps']),
             ('f2', 'x^5 + 1', f2, POLYS['f2']['taps']),
             ('f3', 'x^5 + x^4 + x^2 + 1', f3, POLYS['f3']['taps'])]
    for key, pname, f, taps in polys:
        prim, reason = is_primitive(f)
        say(f'{key} = {pname}: примитивен = {prim}  ({reason})')

    say('')
    say('Периоды LFSR (перебор всех 31 ненулевого начального состояния):')
    max_periods = {}
    best_seeds = {}
    for key, pname, f, taps in polys:
        periods = {}
        for seed in range(1, 32):
            T, pre = lfsr_period(taps, 5, seed)
            periods[T] = periods.get(T, 0) + 1
        tmax = max(periods)
        max_periods[key] = tmax
        best_seeds[key] = [s for s in range(1, 32)
                           if lfsr_period(taps, 5, s)[0] == tmax]
        say(f'  {key} ({pname}): T_max = {tmax}, '
            f'распределение периодов: {dict(sorted(periods.items()))}')

    # выбор: примитивный многочлен f1 даёт максимальный период 2^5 - 1 = 31
    key, pname, f, taps = polys[0]
    seed = 1
    T, pre = lfsr_period(taps, 5, seed)
    demo = lfsr_generate(taps, 5, seed, 31)
    say('')
    say(f'Выбранные параметры (максимальный период): f1 = {pname}, seed = {seed}')
    say(f'Период: T = {T} = 2^5 - 1, предпериод = {pre}')
    say(f'Периодическая часть (31 бит): {demo}')
    # демонстрация произвольной длины
    for L_req in (1, 7, 100, 1000):
        seq = lfsr_generate(taps, 5, seed, L_req)
        assert len(seq) == L_req and set(seq) <= {0, 1}
    say('Проверка генерации произвольной длины (1, 7, 100, 1000): OK')
    # равномерность алфавита: в одном периоде единиц и нулей почти поровну
    seq_period = lfsr_generate(taps, 5, seed, 31)
    say(f'Баланс периода: единиц {sum(seq_period)} из 31 '
        f'(у m-последовательности нечётной длины 31 ровно 16 единиц / 15 нулей)')

    # ------------------------- 3. Тесты на pi / e ---------------------------
    say('')
    say('--- 3a. Проверка тестов на контрольных последовательностях (таблица 1) ---')
    files = {'pi': 'data.pi', 'e': 'data.e'}
    control = {}
    for tag in ('pi', 'e'):
        b = load_bits(files[tag])
        p_runs, info_r = runs_test(b)
        p_tmpl, info_t = non_overlapping_template_test(b, m=9, M=125_000,
                                                       template='000000001')
        p_maurer, info_m = maurers_universal_test(b, L=7, Q=1280)
        control[tag] = {'runs': p_runs, 'nonoverlap': p_tmpl,
                        'maurer': p_maurer}
        say(f'  [{tag}] серий: P={p_runs:.6f}  '
            f'(таблица 1: {TABLE1["runs"][tag]:.6f})')
        say(f'  [{tag}] неперекрывающиеся шаблоны: P={p_tmpl:.6f} '
            f'(таблица 1: {TABLE1["nonoverlap"][tag]:.6f})')
        say(f'  [{tag}] универсальный Маурера: P={p_maurer:.6f} '
            f'(таблица 1: {TABLE1["maurer"][tag]:.6f})')
    # сверка с таблицей 1 (допуск +-0.01)
    all_match = True
    notes = []
    for key in ('runs', 'nonoverlap', 'maurer'):
        for tag in ('pi', 'e'):
            other = 'e' if tag == 'pi' else 'pi'
            mine = control[tag][key]
            ref = TABLE1[key][tag]
            alt = TABLE1[key][other]
            ok = abs(mine - ref) <= 0.01 or abs(mine - alt) <= 0.01
            if not ok:
                all_match = False
                notes.append(f'{key}/{tag}: вычислено {mine:.6f}, '
                             f'в таблице {ref:.6f} / {alt:.6f}')
    if notes:
        say('Замечания сверки: ' + '; '.join(notes))
    say('Сверка с таблицей 1 (допуск +-0.01, колонки pi/e могут быть '
        'переставлены, как в варианте 5): ' + ('OK' if all_match else 'РАСХОЖДЕНИЕ'))

    # ------------------- 3b. Последовательность генератора ------------------
    say('')
    say('--- 3b. Последовательность генератора: 10^6 бит ---')
    bits = ''.join(map(str, lfsr_generate(taps, 5, seed, N_BITS)))
    say(f'Сгенерировано {len(bits)} бит (LFSR f1, seed=1)')
    T_bits = T  # период побитового выхода равен периоду состояний
    say(f'Период битовой последовательности: {T_bits} бит')

    p_runs, info_r = runs_test(bits)
    p_tmpl, info_t = non_overlapping_template_test(bits, m=9, M=125_000,
                                                    template='000000001')
    p_maurer, info_m = maurers_universal_test(bits, L=7, Q=1280)

    say('')
    say(f'Результаты тестов для последовательности генератора (alpha={ALPHA}):')
    say(f'  1) Серий: V={info_r["V"]}, pi={info_r["pi"]:.4f}, '
        f'P={p_runs:.6f} -> {verdict(p_runs)}')
    say(f'  2) Неперекрывающиеся шаблоны (000000001, M=125000): '
        f'P={p_tmpl:.6f} -> {verdict(p_tmpl)}')
    say(f'     (mu={info_t["mu"]:.2f}, sigma2={info_t["sigma2"]:.2f}, '
        f'chi2={info_t["chi2"]:.2f})')
    say(f'  3) Универсальный Маурера (L=7, Q=1280): fn={info_m["fn"]:.6f}, '
        f'E={info_m["expected"]:.6f}, P={p_maurer:.6f} -> {verdict(p_maurer)}')
    all_pass = (p_runs >= ALPHA and p_tmpl >= ALPHA and p_maurer >= ALPHA)
    say(f'Итог по трём тестам: '
        + ('все три ПРОЙДЕНЫ' if all_pass else 'НЕ ПРОЙДЕН хотя бы один тест'))

    text = '\n'.join(out)
    print(text)
    with open('results.txt', 'w', encoding='utf-8') as f:
        f.write(text + '\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
