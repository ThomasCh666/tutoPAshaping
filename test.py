# determine first n prime numbers and print them out  
def get_prime(n):
    prime_list = []
    for i in range(2, n+1):
        for j in range(2, i):
            if i%j == 0:
                break
        else:
            prime_list.append(i)
    return prime_list  
