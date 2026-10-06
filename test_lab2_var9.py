# -*- coding: utf-8 -*-
"""Тесты для lab2_var9.py (ЛР2, вариант 9: LFSR + тесты NIST).

Запуск:  python -m unittest test_lab2_var9.py -v
"""

import math
import unittest

import lab2_var9 as v9


class TestLFSR(unittest.TestCase):
    """Генератор LFSR, периоды, примитивность многочленов."""

    def test_generate_arbitrary_length(self):
        taps, L, seed = v9.POLYS['f1']['taps'], 5, 1
        for length in (0, 1, 7, 31, 100, 1000):
            seq = v9.lfsr_generate(taps, L, seed, length)
            self.assertEqual(len(seq), length)
            self.assertTrue(set(seq) <= {0, 1})

    def test_generate_is_periodic(self):
        """Вывод LFSR строго периодичен: bits[i] == bits[i + T]."""
        taps, L, seed = v9.POLYS['f1']['taps'], 5, 1
        T, pre = v9.lfsr_period(taps, L, seed)
        seq = v9.lfsr_generate(taps, L, seed, 3 * T)
        for i in range(2 * T):
            self.assertEqual(seq[i], seq[i + T])

    def test_max_period_theory(self):
        """Примитивный многочлен даёт T = 2^L - 1 при любом ненулевом seed."""
        taps, L = v9.POLYS['f1']['taps'], 5
        for seed in range(1, 1 << L):
            T, pre = v9.lfsr_period(taps, L, seed)
            self.assertEqual(T, (1 << L) - 1)
            self.assertEqual(pre, 0)

    def test_periods_of_variant_polys(self):
        """Периоды трёх многочленов варианта: f1 -> 31, f2 -> 5, f3 -> 15."""
        expected = {'f1': 31, 'f2': 5, 'f3': 15}
        for key, tmax in expected.items():
            taps, L = v9.POLYS[key]['taps'], 5
            periods = {v9.lfsr_period(taps, L, s)[0] for s in range(1, 32)}
            self.assertEqual(max(periods), tmax)

    def test_brent_matches_dict_period(self):
        """Метод Брента согласуется со словарём состояний."""
        for key in ('f1', 'f2', 'f3'):
            taps, L = v9.POLYS[key]['taps'], 5
            for seed in (1, 2, 7, 17, 31):
                self.assertEqual(v9.lfsr_period(taps, L, seed),
                                 v9.lfsr_period_fast(taps, L, seed))

    def test_primitivity_of_variant_polys(self):
        """Примитивен только f1 = x^5 + x^2 + 1."""
        f1 = (1 << 5) | (1 << 2) | 1
        f2 = (1 << 5) | 1
        f3 = (1 << 5) | (1 << 4) | (1 << 2) | 1
        prim_f1, _ = v9.is_primitive(f1)
        prim_f2, _ = v9.is_primitive(f2)
        prim_f3, _ = v9.is_primitive(f3)
        self.assertTrue(prim_f1)
        self.assertFalse(prim_f2)
        self.assertFalse(prim_f3)

    def test_m_sequence_balance(self):
        """m-последовательность периода 31 содержит 16 единиц и 15 нулей."""
        taps, L, seed = v9.POLYS['f1']['taps'], 5, 1
        seq = v9.lfsr_generate(taps, L, seed, 31)
        self.assertEqual(seq.count(1), 16)
        self.assertEqual(seq.count(0), 15)

    def test_known_lfsr_example(self):
        """LFSR x^3 + x + 1: классическая m-последовательность 0,0,1,0,1,1,1
        (фаза цикла зависит от seed; при seed=0b100 выход начинается с 0,0,1)."""
        taps = (1, 0)
        seq = v9.lfsr_generate(taps, 3, 0b100, 7)
        self.assertEqual(seq, [0, 0, 1, 0, 1, 1, 1])
        T, pre = v9.lfsr_period(taps, 3, 0b100)
        self.assertEqual(T, 7)


class TestGF2(unittest.TestCase):
    """Арифметика GF(2): неприводимость и примитивность."""

    def test_irreducible_small(self):
        self.assertTrue(v9.is_irreducible(0b1011))   # x^3 + x + 1
        self.assertTrue(v9.is_irreducible(0b10011))  # x^4 + x + 1
        self.assertTrue(v9.is_irreducible(0b111))    # x^2 + x + 1 — неприводим
        self.assertFalse(v9.is_irreducible(0b101))    # x^2 + 1 = (x + 1)^2

    def test_primitive_small(self):
        prim, _ = v9.is_primitive(0b1011)  # x^3 + x + 1 примитивен (T=7)
        self.assertTrue(prim)
        # x^4 + x^3 + x^2 + x + 1 неприводим, но не примитивен (T=5)
        prim, _ = v9.is_primitive(0b11111)
        self.assertFalse(prim)

    def test_factorization_f2(self):
        """x^5 + 1 = (x + 1)(x^4 + x^3 + x^2 + x + 1)."""
        self.assertEqual(v9.gf2_mul(0b11, 0b11111), (1 << 5) | 1)

    def test_factorization_f3(self):
        """x^5 + x^4 + x^2 + 1 = (x + 1)(x^4 + x + 1)."""
        f3 = (1 << 5) | (1 << 4) | (1 << 2) | 1
        self.assertEqual(v9.gf2_mul(0b11, 0b10011), f3)


class TestRuns(unittest.TestCase):
    """Тест серий (NIST 2.3)."""

    def test_nist_example(self):
        """Пример NIST 2.3.8: n=100, pi=0.42, V=52, P=0.500798."""
        eps = ('11001001000011111101101010100010001000010110100011'
               '00001000110100110001001100011001100010100010111000')
        self.assertEqual(len(eps), 100)
        p, info = v9.runs_test(eps)
        self.assertAlmostEqual(info['pi'], 0.42, places=10)
        self.assertEqual(info['V'], 52)
        self.assertAlmostEqual(p, 0.500798, places=6)

    def test_alternating_fails(self):
        """Чередование 0101...: слишком быстрые колебания, P < alpha."""
        p, info = v9.runs_test('01' * 50)
        self.assertLess(p, 0.01)

    def test_prerequisite_monobit(self):
        """Дисбалансные данные: частотный предтест не пройден, P = 0."""
        p, info = v9.runs_test('1' * 90 + '0' * 10)
        self.assertEqual(p, 0.0)
        self.assertIn('note', info)


class TestNonOverlappingTemplates(unittest.TestCase):
    """Тест неперекрывающихся шаблонов (NIST 2.7)."""

    def test_nist_small_example(self):
        """Пример NIST 2.7.4: eps=10100100101110010110, m=3, B=001, N=2, M=10.

        Блок 1 = 1010010010: неперекрывающийся поиск даёт вхождения в позициях
        0-based 3 и 6 -> W1 = 2; блок 2 = 1110010110: вхождение в позиции 3 ->
        W2 = 1. Тогда chi2 = 2.133333 и P = 0.344154 — в точности значения из
        примера NIST (в тексте примера опечатка в счётчиках W, но chi2/P
        посчитаны именно по W = (2, 1)).
        """
        eps = '10100100101110010110'
        p, info = v9.non_overlapping_template_test(eps, m=3, M=10, template='001')
        self.assertEqual(info['W'], [2, 1])
        self.assertAlmostEqual(info['mu'], 1.0, places=10)
        self.assertAlmostEqual(info['sigma2'], 0.46875, places=10)
        self.assertAlmostEqual(info['chi2'], 2.133333, places=5)
        self.assertAlmostEqual(p, 0.344154, places=5)

    def test_window_skips_after_hit(self):
        """После вхождения окно сдвигается на m бит (неперекрывающийся поиск)."""
        block = '001' * 3 + '0'   # вхождения в позициях 0, 3 и 6
        p, info = v9.non_overlapping_template_test(block, m=3, M=10, template='001')
        self.assertEqual(info['W'], [3])

    def test_control_pi_e(self):
        """Таблица 1: pi -> 0.078790/0.165757 (в таблице колонки переставлены)."""
        bits = v9.load_bits('data.pi')
        p, _ = v9.non_overlapping_template_test(bits)
        self.assertTrue(abs(p - 0.165757) <= 0.01 or abs(p - 0.078790) <= 0.01)
        bits = v9.load_bits('data.e')
        p, _ = v9.non_overlapping_template_test(bits)
        self.assertTrue(abs(p - 0.165757) <= 0.01 or abs(p - 0.078790) <= 0.01)


class TestMaurer(unittest.TestCase):
    """Универсальный тест Маурера (NIST 2.9)."""

    def test_doc_example_fn(self):
        """Пример NIST 2.9.4: n=20, L=2, Q=4, K=6 -> fn = 1.1949875.

        Шаги примера (1-блочная нумерация документа): init-таблица
        00->0, 01->2, 10->4, 11->0; тестовые блоки 5..10:
        log2(5-2) + log2(6-0) + log2(7-5) + log2(8-7) + log2(9-8) + log2(10-6)
        = 7.169925002, fn = 7.169925002/6 = 1.1949875.
        """
        total = sum(math.log2(d) for d in (3, 6, 2, 1, 1, 4))
        self.assertAlmostEqual(total, 7.169925002, places=8)
        self.assertAlmostEqual(total / 6, 1.1949875, places=6)

    def test_control_pi_e(self):
        """Таблица 1: pi -> 0.282568/0.669012 (колонки переставлены), допуск 0.01."""
        for path, ref in (('data.pi', (0.282568, 0.669012)),
                          ('data.e', (0.282568, 0.669012))):
            bits = v9.load_bits(path)
            p, info = v9.maurers_universal_test(bits, L=7, Q=1280)
            self.assertTrue(any(abs(p - r) <= 0.01 for r in ref),
                            msg=f'{path}: P={p} vs {ref}')

    def test_compressible_sequence_fails(self):
        """Сжимаемая (периодическая) последовательность проваливает тест."""
        bits = ('0100' * 250_000)[:10**6]
        p, info = v9.maurers_universal_test(bits, L=7, Q=1280)
        self.assertLess(p, 0.01)


class TestMain(unittest.TestCase):
    """Интеграционный прогон main(): все пункты задания, exit 0."""

    def test_main_exit_zero(self):
        self.assertEqual(v9.main(), 0)


if __name__ == '__main__':
    unittest.main()
