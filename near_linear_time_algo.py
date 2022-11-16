from math import log2


k = 1
rho_p_i = 0
rho_p_i_prime = 0
rho_p_d = 0
rho_p_d_prime = 0
rho_p_s = 0


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
    f = {0: 0, 1: 1}

    J = 2 * (max(1, round(3/2 * (rho_p_i * rho_p_i_prime +
             (rho_p_d + 1/(k * log2(len(s1)))) * (rho_p_d_prime + 1))) + 1) * k)
    I = round(len(s1)/(k*log2(len(s1)))) - 1

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
                f[round(i*k*log2(len(s1)) + 1)] = f.get(round((i - 1)*k *
                                                              log2(len(s2))), 0) + round((j + k)*log2(len(s2)))

    return f


app = approx_align(['C', 'A', 'G', 'G', 'G', 'C', 'A', 'G', 'G', 'G', 'C', 'A', 'G', 'G',
                    'G', 'C', 'A', 'G', 'G', 'G', 'C', 'A', 'G', 'G', 'G', 'C', 'A', 'G', 'G', 'G'], ['C', 'A', 'G', 'G', 'G', 'G', 'A', 'G', 'G', 'G', 'A', 'G', 'G', 'G', 'G', 'A', 'G', ])

print(app)
