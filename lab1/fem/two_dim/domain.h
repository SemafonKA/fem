#pragma once
#include <vector>



namespace fem::two_dim {

    struct Subdomain {
        size_t materialNum{ 0 }; /// Number of material
        size_t xBeginNum{ 0 };   /// Number of left X coordinate line 
        size_t xEndNum{ 0 };     /// Number of right X coordinate line
        size_t yBeginNum{ 0 };   /// Number of bottom Y coordinate line 
        size_t yEndNum{ 0 };     /// Number of top Y coordinate line 
    };

    /**
     * @brief Description of domain of selected field
     */
    struct Domain {
        size_t Kx{}; /// count of X coordinate lines
        size_t Ky{}; /// count of Y coordinate lines

        std::vector<double> X{}; /// Vector of X coordinate lines of Kx*Ky size
        std::vector<double> Y{}; /// Vector of Y coordinate lines of Kx*Ky size

        std::vector<Subdomain> subdomains{}; /// Subdomains of current domain

        std::vector<size_t> nx{}; /// Count of subdivides for X coordinate lines of (Kx*Ky) - 1 size
        std::vector<double> cx{}; /// Sparse coefficients for X coordinate lines of (Kx*Ky) - 1 size
        std::vector<size_t> ny{}; /// Count of subdivides for Y coordinate lines of (Kx*Ky) - 1 size
        std::vector<double> cy{}; /// Sparse coefficients for Y coordinate lines of (Kx*Ky) - 1 size

        size_t splitX{ 0 }; /// Count of additional subdivides for X coordinate lines (0 is no additional subdivides)
        size_t splitY{ 0 }; /// Count of additional subdivides for Y coordinate lines (0 is no additional subdivides)
    };

}
