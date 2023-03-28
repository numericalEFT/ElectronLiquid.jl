# Exact diagonalization of Fermi Hubbard model
using LinearAlgebra

function fermi_hubbard_model(N, t, U)
    # N: number of lattice sites
    # t: hopping parameter
    # U: on-site interaction

    # Initialize Hamiltonian matrix
    # The Hamiltonian matrix is a 2^N x 2^N matrix, where N is the number of lattice sites.
    # This is because there are 2^N possible basis states for the system.
    H = zeros(2^N, 2^N)

    # Loop over all basis states
    # The loop starts from 0 because the basis states are represented as binary numbers,
    # and the first basis state is all zeros.
    for i in 0:(2^N - 1)
        # Loop over all lattice sites
        for j in 1:N
            # Calculate hopping term
            # The hopping term is calculated by XORing the current basis state (i) with
            # the binary representation of the lattice site (j-1) and its neighbor.
            # This is done for both up and down spins (j-1 + N for down spins).
            # The hopping term is subtracted from the diagonal element of the Hamiltonian matrix.
            if j < N
                target_up = xor(i, (1 << (j-1))) + xor(i, (1 << j))
                target_down = xor(i, (1 << (j-1 + N))) + xor(i, (1 << (j + N)))
                H[i+1, target_up+1] -= t
                H[i+1, target_down+1] -= t
            else
                target_up = xor(i, (1 << (j-1))) + xor(i, (1 << (j-1 - N)))
                target_down = xor(i, (1 << (j-1 + N))) + xor(i, (1 << (j-1)))
                H[i+1, target_up+1] -= t
                H[i+1, target_down+1] -= t
            end

            # Calculate on-site interaction term
            # The on-site interaction term is calculated by checking if both up and down spins
            # are occupied at the same lattice site (j-1 for up spins and j-1 + N for downspins).
            # If both spins are occupied, the on-site interaction term (U) is added to the
            # diagonal element of the Hamiltonian matrix.
            if (i & (1 << (j-1)) != 0) && (i & (1 << (j-1 + N)) != 0)
                H[i+1, i+1] += U
            end
        end
    end    # After looping through all basis states and lattice sites, the Hamiltonian matrix (H)
    # is fully constructed with the hopping and on-site interaction terms.
    # The function then returns the Hamiltonian matrix.
    return H
end

# Example usage
N = 4  # Number of lattice sites
t = 1.0  # Hopping parameter
U = 2.0  # On-site interaction
H = fermi_hubbard_model(N, t, U)

# Diagonalize the Hamiltonian
eigenvalues, eigenvectors = eigen(H)
