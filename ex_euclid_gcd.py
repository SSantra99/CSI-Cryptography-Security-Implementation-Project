def ex_euclid_gcd(a, b):
    old_s, s = 1, 0
    old_t, t = 0, 1
    table=[]
    while b != 0:
        quotient = a//b
        a,b=b,a%b
        old_s, s = s, old_s - quotient * s
        old_t, t = t, old_t - quotient * t
        table.append([a,quotient,old_s,old_t])

    print("Final coefficient of the given numbers to get their gcd:", (old_s, old_t))
    print("Greatest common divisor:", a)
    print(table)


# Example usage
a = 306
b = 189
ex_euclid_gcd(a, b)
