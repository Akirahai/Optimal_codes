from itertools import product
from collections import defaultdict
from tabulate import tabulate
from tqdm import tqdm  # For progress 
import math

# Optimal Codes

class Vertex:
    """Class representing a vertex in the q-ary graph."""
    def __init__(self, sequence):
        self.sequence = sequence  # Store the sequence (e.g., '001')
        self.Bx = self.generate_Bx()  # Store the B(x) set for the vertex

    def generate_Bx(self):
        """Generate the set B(x) by removing each position from the sequence."""
        Bx = set()
        for i in range(len(self.sequence)):
            Bx.add(self.sequence[:i] + self.sequence[i + 1:])
        return Bx

    def __str__(self):
        return f"Vertex({self.sequence})"

    def __repr__(self):
        return self.__str__()

class QAryGraph:
    """Graph of q-ary sequences with length n as vertices."""
    def __init__(self, q, n):
        self.q = q  # q- ary sequences
        self.n = n  # Length of the sequences
        self.vertices = self.generate_vertices()  # List of Vertex objects
        self.adj_list = defaultdict(list)  # Adjacency list
        self.adj_matrix = [[0] * len(self.vertices) for _ in range(len(self.vertices))]
        self.vertex_indices = {v.sequence: i for i, v in enumerate(self.vertices)}  # Map sequence to index

    def generate_vertices(self):
        """Generate all q-ary sequences of length n as Vertex objects."""
        sequences = [''.join(map(str, seq)) for seq in product(range(self.q), repeat=self.n)]
        return [Vertex(sequence) for sequence in sequences]

    def are_connected(self, v1, v2):
        """Check if two vertices are connected (B(x) ∩ B(y) = ∅)."""
        return v1.Bx.isdisjoint(v2.Bx)

    def build_graph(self):
        """Build the graph by constructing the adjacency list and matrix."""
        for i, v1 in enumerate(self.vertices):
            for j, v2 in enumerate(self.vertices):
                if i < j and self.are_connected(v1, v2):
                    # Add to adjacency list
                    self.adj_list[v1.sequence].append(v2.sequence)
                    self.adj_list[v2.sequence].append(v1.sequence)
                    # Update adjacency matrix
                    self.adj_matrix[i][j] = 1
                    self.adj_matrix[j][i] = 1

    def get_adjacency_list(self):
        """Return the adjacency list."""
        return dict(self.adj_list)

    def get_adjacency_matrix(self):
        """Return the adjacency matrix."""
        return self.adj_matrix

    def bron_kerbosch(self, R, P, X, cliques,progress, total_calls, call_counter):
        """Bron-Kerbosch algorithm to find all maximal cliques."""
        progress.update(1)
        call_counter[0] += 1  # Track the number of recursive calls

        # Update progress every 25%
        if call_counter[0] >= total_calls // 4 * progress.n + 1:
            progress.update(total_calls // 4)
            
        if not P and not X:
            cliques.append(R)
            return
        
        for v in list(P):
            neighbors = set(self.adj_list[v])
            self.bron_kerbosch(R | {v}, P & neighbors, X & neighbors, cliques,progress, total_calls, call_counter)
            P.remove(v)
            X.add(v)

    def max_clique(self):
        """Find all maximum cliques and their size."""
        cliques = []
        P = set(v.sequence for v in self.vertices)
        
        total_vertices = len(self.vertices)
        total_calls = math.ceil(3 ** (total_vertices / 3))
        
        with tqdm(total=total_calls, desc=f"Finding optimal codes for {self.q}- ary sequences with length {self.n}", unit="call") as progress:
            call_counter = [0]
            self.bron_kerbosch(set(), P, set(), cliques, progress, total_calls, call_counter)

        # Find the size of the maximum cliques
        max_size = max(len(clique) for clique in cliques)

        # Collect all cliques with the maximum size
        max_cliques = [clique for clique in cliques if len(clique) == max_size]

        return max_size, max_cliques
    
    
    
# Formula to find VT_0(n) module n + 1

from math import gcd, isqrt

class VT_Formula:
    def __init__(self, a, n):
        """
        Initialize the VT_Formular class with parameters a and n.
        """
        self.a = a
        self.n = n

    def phi(self, n):
        """
        Euler's Totient Function φ(n)
        Counts the number of integers from 1 to n that are coprime with n.
        """
        result = n
        p = 2
        while p * p <= n:
            if n % p == 0:
                # Subtract multiples of p
                while n % p == 0:
                    n //= p
                result -= result // p
            p += 1
        if n > 1:
            result -= result // n
        return result

    def mobius(self, n):
        """
        Möbius Function μ(n)
        Returns:
        - 0 if n has squared prime factors,
        - (-1)^k where k is the number of distinct primes dividing n, otherwise.
        """
        if n == 1:
            return 1
        primes_count = 0
        for p in range(2, isqrt(n) + 1):
            if n % p == 0:
                if n % (p * p) == 0:
                    return 0  # n has squared prime factor
                primes_count += 1
                n //= p
            while n % p == 0:
                n //= p
        if n > 1:
            primes_count += 1
        return -1 if primes_count % 2 else 1

    def c(self, d,a):
        """
        Helper function c(d, a) used in V_Ta_n formula.
        """
        phi_d = self.phi(d)
        mu_term = self.mobius(d // gcd(d, a))
        phi_term = self.phi(d // gcd(d, a))
        return phi_d * (mu_term // phi_term)

    def V_Ta_n(self):
        """
        Main function to calculate the |V_Ta(n)| value based on the provided formula.
        """
        total = 0
        n_plus_1 = self.n + 1

        # Find all divisors d of n+1, such that d is odd
        for d in range(1, n_plus_1 + 1):
            if n_plus_1 % d == 0 and d % 2 == 1:  # d divides n+1 and d is odd
                total += self.c(d,self.a) * (2 ** ((self.n + 1) // d))

        return total // (2 * self.n + 2)
    
    
    
    
    
    
# Formula to find VT(a,q,n) module q*n

class VT:
    def __init__(self, a, q, n, module_value):
        """
        Initialize the VT class with a, q, and n.
        a: the target mod value
        q: the upper bound for the variables x_i (x_i ∈ {0, ..., q-1})
        n: the number of variables in the equation
        """
        self.a = a  # target value for S mod q*n
        self.q = q  # upper bound for x_i values (x_i ∈ {0, ..., q-1})
        self.n = n  # number of variables in the equation
        self.module = module_value  # modulus value q*n

    def find_solutions(self, a):
        """
        Finds the number of solutions {x_1, x_2, ..., x_n} 
        such that the equation:
        x_1 + 2*x_2 + ... + n*x_n = S where S ≡ a (mod q*n) holds.
        """
        mod_value = self.module  # modulus value q*n
        count = 0

        # Generate all possible combinations of {x_1, ..., x_n} where x_i ∈ {0, ..., q-1}
        for x in product(range(self.q), repeat=self.n):
            # Calculate S = x_1 + 2*x_2 + ... + n*x_n
            S = sum((i + 1) * x[i] for i in range(self.n))
            a = a % mod_value    
            # Check if S ≡ a (mod q*n)
            if S % mod_value == a:
                count += 1

        return count

    def w_a(self, a):
        """
        Returns the value of w_a(q, n) which is the number of solutions for a.
        """
        return self.find_solutions(a)
    
    
    def print_all_solutions(self):
        """
        Prints the number of solutions for each value of a from 0 to q*n - 1.
        """
        mod_value = self.module # modulus value q*n

        total_solutions = 0
        for a in range(mod_value):
            num_solutions = self.find_solutions(a)
            print(f"a = {a}: Number of solutions = {num_solutions}")
            total_solutions += num_solutions
        
        print(f"Total number of solutions = {total_solutions}")
        print(f'Total number of solutions checking again: {self.q ** self.n}')
        
    def compute_fz(self, z):
        """
        Computes f(z) = sum of w_a(q, n) * z^a for a in range(0, self.module).
        z: the variable in the function f(z)
        """
        mod_value = self.module  # modulus value q*n
        fz = 0

        # Calculate f(z) = Sum of w_a(q, n) * z^a for a from 0 to q*n - 1
        for a in range(mod_value + 1):
            w_a_value = self.w_a(a)  # Find w_a(q, n)
            fz += w_a_value * (z ** a)  # Add w_a(q, n) * z^a to the sum

        return fz
            