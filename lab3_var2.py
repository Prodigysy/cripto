#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Лабораторная работа № 3. Генерация простого числа.
Вариант 2 — вероятностный тест простоты Леманна.

Схема генерации (по методичке):
  1. Генерируется случайное p-битовое число n;
  2. старший бит ставится в 1 (гарантия длины p бит),
     младший бит ставится в 1 (гарантия нечётности);
  3. n проверяется на делимость на маленькие простые числа
     (все простые меньше 2000 — наиболее эффективный вариант по методичке);
  4. выполняется 5 раундов теста Леманна со случайными a.
     Если n не проходит хотя бы один раунд — берётся новый кандидат.

Тест Леманна (раздел 3.6 методички):
  1. Выбрать случайное a < n;
  2. вычислить t = a^((n-1)/2) mod n;
  3. если t != ±1 (mod n), то n — составное;
  4. иначе ответ «неизвестно», вероятность ошибки не больше 1/2.
  После k повторений вероятность пропустить составное число <= 1/2^k.

Запуск:
    python lab3_var2.py            # p = 256 по умолчанию
    python lab3_var2.py 512        # другая длина в битах

Зависимости: только стандартная библиотека Python 3.
"""

import random
import sys

# Предел малых простых для пробных делений (по методичке — «меньше 2000»)
SMALL_PRIME_LIMIT = 2000
# Число раундов теста Леманна (по методичке — 5)
LEHMANN_ROUNDS = 5
# Верхняя граница для случайного a: методичка рекомендует небольшие
# значения a для ускорения вычислений
A_BOUND = 10 ** 6


def sieve(limit):
    """Решето Эратосфена: все простые <= limit."""
    is_prime = [True] * (limit + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(limit ** 0.5) + 1):
        if is_prime[i]:
            for j in range(i * i, limit + 1, i):
                is_prime[j] = False
    return [i for i, ok in enumerate(is_prime) if ok]


# Все простые < 2000 для быстрой предварительной отбраковки
SMALL_PRIMES = sieve(SMALL_PRIME_LIMIT - 1)


def generate_candidate(p):
    """Случайное p-битовое нечётное число: старший и младший биты равны 1."""
    n = random.getrandbits(p)
    n |= 1 << (p - 1)   # старший бит — гарантия длины p бит
    n |= 1              # младший бит — гарантия нечётности
    return n


def small_prime_check(n):
    """Пробные деления на простые < 2000. Возвращает делитель или None."""
    for p in SMALL_PRIMES:
        if n % p == 0:
            return p
    return None


def lehmann_test(n, a):
    """Один раунд теста Леманна с основанием a.

    Возвращает True, если a^((n-1)/2) = ±1 (mod n) — «неизвестно»,
    False, если найдено основание, выдающее составность n.
    """
    t = pow(a, (n - 1) // 2, n)
    return t == 1 or t == n - 1


def generate_prime(p, verbose=True):
    """Генерация простого p-битового числа по схеме методички."""
    candidates = 0          # всего кандидатов
    rejected_small = 0      # отсеяно пробными делениями
    rejected_lehmann = 0    # отсеяно тестом Леманна

    while True:
        candidates += 1
        n = generate_candidate(p)

        # Шаг 3: делимость на маленькие простые
        d = small_prime_check(n)
        if d is not None:
            rejected_small += 1
            if verbose:
                print(f"  кандидат #{candidates}: {n}\n"
                      f"    отбракован пробным делением: n % {d} == 0")
            continue

        # Шаг 4: 5 раундов теста Леманна
        if verbose:
            print(f"  кандидат #{candidates}: {n}\n"
                  f"    пробные деления пройдены, тест Леманна:")
        prime = True
        for round_num in range(1, LEHMANN_ROUNDS + 1):
            a = random.randint(2, min(n - 2, A_BOUND))
            ok = lehmann_test(n, a)
            if verbose:
                t = pow(a, (n - 1) // 2, n)
                print(f"    раунд {round_num}: a = {a}, "
                      f"a^((n-1)/2) mod n = {t} -> "
                      f"{'±1, неизвестно' if ok else 'не ±1, составное'}")
            if not ok:
                prime = False
                rejected_lehmann += 1
                break
        if prime:
            return n, candidates, rejected_small, rejected_lehmann
        if verbose:
            print("    n отвергнуто, генерируем новый кандидат")


def main():
    p = int(sys.argv[1]) if len(sys.argv) > 1 else 256
    print(f"ЛР №3, вариант 2: генерация простого {p}-битового числа "
          f"(тест Леманна, {LEHMANN_ROUNDS} раундов)\n")

    print(f"[1] Решето Эратосфена: {len(SMALL_PRIMES)} простых < "
          f"{SMALL_PRIME_LIMIT} для пробных делений\n")

    print(f"[2] Генерация простого {p}-битового числа:")
    n, candidates, rejected_small, rejected_lehmann = generate_prime(p)

    print(f"\n[3] Итог:")
    print(f"    кандидатов рассмотрено:      {candidates}")
    print(f"    отбраковано пробными делениями: {rejected_small}")
    print(f"    отбраковано тестом Леманна:  {rejected_lehmann}")
    print(f"    ПРОСТОЕ ЧИСЛО ({p} бит):\n{n}")
    print(f"\n    Вероятность ошибки после {LEHMANN_ROUNDS} раундов "
          f"теста Леманна <= 1/2^{LEHMANN_ROUNDS} = {1.0 / 2 ** LEHMANN_ROUNDS}")


if __name__ == "__main__":
    main()
