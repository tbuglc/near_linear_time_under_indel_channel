from math import log2, ceil

influenza = "AGCAAAAGCAGGAGTTCAAAATGAATCCAAATCAGAAAATAATAACCATTGGGTCAATCTGTATGGGGATCGGAATAATCAGCCTAATGCTCTCACTTGGAATGCAGGACCTTTTTTCTGACTCAGGGCGCCTTGTTGAATGACAAACACTCAAATGGAACCGTTAAAGACAGAAGCCCTTATAGAACCTTGATGAGCTGTCCTGTTGGTGAAGCTCCTTCTCCCTACAATTCAAGGTTCGGTGT"

monkeypox = "CTCTTTCTCTCTTCGATGGGTCTCACAAAAATATTAAACCAGAATATATTGTTGGACGTTATCGTTTACGAAATAGTTGAGACATCAGAAAGAGGTTTAATATTTTTGTGAGACCATCGAAGAGAGAAAGAGAATAAAAATATTCGATTAACCCAACTCATCCATTTTCAGATGAATAGAGTTATCGATTCAGACACATGCTTTGAGTTTTGTTGAATCGATGAGTGAAGTATCATCGGTTGCACCTTCAGATC"
# TODO: PARAMS CHANGES
k = 10
rho_p_i = .028
rho_p_i_prime = .028
rho_p_d = .028
rho_p_d_prime = .028
rho_p_s = 0.34


def get_approximate_distance(apprx, s1, s2):
    s_p_1 = ''
    s_p_2 = ''

    for key in apprx:
        s_p_1 += s1[slice(key - 1, key + ceil(k*log2(len(s1))) - 1)]
        s_p_2 += s2[slice(apprx[key] - 1, apprx[key] +
                          ceil(k*log2(len(s2))) - 1)]

    return s1, s2


def ED(str1, str2):
    m = len(str1)
    n = len(str2)

    dp = [[0 for x in range(n + 1)] for x in range(m + 1)]
    for i in range(m + 1):
        for j in range(n + 1):
            if i == 0:
                dp[i][j] = j
            elif j == 0:
                dp[i][j] = i
            elif str1[i-1] == str2[j-1]:
                dp[i][j] = dp[i-1][j-1]
            else:
                dp[i][j] = 1 + min(dp[i][j-1],
                                   dp[i-1][j],
                                   dp[i-1][j-1])

    return dp[m][n]


def approx_align(s1, s2):
    f = {1: 1}

    J = ceil(2 * (max(1, (3/2 * (rho_p_i * rho_p_i_prime +
             (rho_p_d + 1/(k * log2(len(s1)))) * (rho_p_d_prime + 1))) + 1) * k))

    I = ceil(len(s1)/(k*log2(len(s1)))) - 1

    for i in range(1, I):
        minEd = 100000
        for j in range(J, round(J/2), -1):
            s1_1 = s1[
                slice(
                    round(i*k*log2(len(s1)) + 1),

                    round(
                        (i + 1)*k*log2(len(s1))
                    )
                )
            ]

            s2_2 = s2[
                slice(
                    f.get(round((i - 1)*k*log2(len(s2) + 1)), 0) +
                    round((j + k)*log2(len(s2))),
                    f.get(round((i - 1)*k*log2(len(s2) + 1)), 0) +
                    round((j+2*k) * log2(len(s2)) - 1)
                )
            ]

            ed = ED(s1_1, s2_2)

            if (ed <= minEd):
                minEd = ed
                # f[round(i*k*log2(len(s1))) + 1] =
                f[round(i*k*log2(len(s1)))] = f.get(round((i - 1)*k *
                                                          log2(len(s2))), 0) + round((j + k)*log2(len(s2)))

    return f


app = approx_align(influenza, monkeypox)
print(app)

# ap_s1, ap_s2 = get_approximate_distance(app, influenza, monkeypox)

# print(ap_s1)
# print('============')
# print(ap_s2)
