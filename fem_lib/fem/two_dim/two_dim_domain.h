#pragma once
#include <vector>
#include <string>

namespace fem::two_dim {

    struct Subdomain {
        size_t materialNum{ 0 }; /// Number of material
        size_t xBeginNum{ 0 };   /// Number of left X coordinate line (starts with 1)
        size_t xEndNum{ 0 };     /// Number of right X coordinate line (starts with 1)
        size_t yBeginNum{ 0 };   /// Number of bottom Y coordinate line (starts with 1) 
        size_t yEndNum{ 0 };     /// Number of top Y coordinate line (starts with 1) 
    };

    struct S1Edge {
        size_t xBeginNum{ 0 }; /// Number of left X coordinate line (starts with 1)
        size_t xEndNum{ 0 };   /// Number of right X coordinate line (starts with 1)
        size_t yBeginNum{ 0 }; /// Number of bottom Y coordinate line (starts with 1) 
        size_t yEndNum{ 0 };   /// Number of top Y coordinate line (starts with 1) 
        size_t funcNum{ 0 };   /// Number of function for this edge condition (starts with 1)
    };

    struct S2Edge {
        size_t xBeginNum{ 0 }; /// Number of left X coordinate line (starts with 1)
        size_t xEndNum{ 0 };   /// Number of right X coordinate line (starts with 1)
        size_t yBeginNum{ 0 }; /// Number of bottom Y coordinate line (starts with 1) 
        size_t yEndNum{ 0 };   /// Number of top Y coordinate line (starts with 1) 
        size_t funcNum{ 0 };   /// Number of function for this edge condition (starts with 1)
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

        std::vector<S1Edge> s1_edges{}; /// S1 Edge conditions
        std::vector<S2Edge> s2_edges{}; /// S2 Edge conditions

        std::vector<size_t> nx{}; /// Count of subdivides for X coordinate lines of (Kx*Ky) - 1 size
        std::vector<double> cx{}; /// Sparse coefficients for X coordinate lines of (Kx*Ky) - 1 size
        std::vector<size_t> ny{}; /// Count of subdivides for Y coordinate lines of (Kx*Ky) - 1 size
        std::vector<double> cy{}; /// Sparse coefficients for Y coordinate lines of (Kx*Ky) - 1 size

        size_t splitX{ 0 }; /// Count of additional subdivides for X coordinate lines (0 is no additional subdivides)
        size_t splitY{ 0 }; /// Count of additional subdivides for Y coordinate lines (0 is no additional subdivides)

        /**
         * @brief Read domain structure like described in example from file (ignoring comments)
         * @param filepath - path to file with domain describing
         * @return Domain structure
         * @throws std::runtime_error - when the file cannot been opened or the data in the file contains errors
         */
        static auto readFromFile(const std::string& filepath) -> Domain;

        /**
         * @brief Make a dump of current domain includes additional comments
         * @return - string representation of dump file
         */
        [[nodiscard]]
        auto dump() -> std::string;

        /**
         * @brief Make a dump of current domain without additional comments
         * @return - string representation of dump file
         */
        [[nodiscard]]
        auto dumpNoComments() -> std::string;
    };

}
