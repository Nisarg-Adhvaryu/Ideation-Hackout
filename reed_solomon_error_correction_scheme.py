import galois

GF16 = galois.GF(2**4)

# generator polynomial 
g = galois.Poly([1, 1, 0, 1], field=GF16)  # g(x) = x^3 + x + 1

def rs_encode(message, n, k):
    # Create the message polynomial
    m = galois.Poly(message, field=GF16)
    
    # Multiply the message polynomial by x^(n-k) 
    m_shifted = m * galois.Poly([1] + [0] * (n - k), field=GF16)
    
    # Compute the remainder
    remainder = m_shifted % g
    
    #codeword 
    codeword = list(m.coeffs) + list(remainder.coeffs)
    
    return codeword


message = [1, 0, 1]  
n = 7  
k = 3  
codeword = rs_encode(message, n, k)
print("Encoded Codeword:", codeword)

#Syndrome Calculation

def rs_syndromes(received, n, t):
    syndromes = [galois.Poly(received, field=GF16)(GF16(2)**i) for i in range(1, 2*t+1)]
    return syndromes

#Berlekamp-Massey Algorithm

def berlekamp_massey(syndromes, t):
    sigma = galois.Poly([1], field=GF16)  
    b = galois.Poly([1], field=GF16)
    
    L = 0  
    m = 1  
    for n in range(len(syndromes)):
        
        discrepancy = syndromes[n]
        for i in range(1, L+1):
            discrepancy += sigma.coeffs[i] * syndromes[n-i]  
        
        if discrepancy == 0:
            m += 1
        else:
            temp = galois.Poly(sigma.coeffs, field=GF16) 
            
            
            leading_coeff = sigma.coeffs[0] 
            
            
            sigma += (discrepancy * b * galois.Poly([1] + [0]*m, field=GF16)) // leading_coeff
            if 2*L <= n:
                L = n + 1 - L
                b = temp // discrepancy
                m = 1
            else:
                m += 1
    
    return sigma

# Chien Search
def chien_search(sigma, n):
    error_positions = []
    for i in range(n):
        if sigma(GF16(2)**i) == 0:
            error_positions.append(n - 1 - i)
    return error_positions

# Forney's Algorithm
def forney_algorithm(syndromes, error_positions, sigma):
    n = len(syndromes)
    
    omega = galois.Poly(syndromes, field=GF16) * sigma
    
    error_values = []
    for pos in error_positions:
        x_inv = GF16(2)**(pos * (-1))
        error_value = omega(x_inv) // sigma.derivative()(x_inv)
        error_values.append(error_value)
    
    return error_values

#error correction

def correct_errors(received, error_positions, error_values):
    corrected = received.copy()
    for i, pos in enumerate(error_positions):
        corrected[pos] -= GF16(error_values[i]) 
    return corrected

#Example run
received = GF16([1, 2, 3, 4, 5, 6, 7])  # Received codeword (with potential errors)
n = 7  
k = 3  
t = 2  

syndromes = rs_syndromes(received, n, t)
sigma = berlekamp_massey(syndromes, t)
error_positions = chien_search(sigma, n)
error_values = forney_algorithm(syndromes, error_positions, sigma)
corrected_codeword = correct_errors(received, error_positions, error_values)

print("Corrected Codeword:", corrected_codeword)
