#include "two_dim_domain.h"

#include <sstream>
#include <string>
#include <format>
#include <algorithm>
#include <vector>

#include "two_dim_domain_reader.h"

using std::ostringstream;
using std::string;
using std::format;
using std::vector;

static
size_t findMaxSize(const vector <double>& vec) {
    size_t max = 0;
    for (const auto& el : vec) {
        max = std::max(max, format("{}", el).size());
    }
    return max;
}

namespace fem::two_dim {

    /**
     * @brief Read domain structure like described in example from file (ignoring comments)
     * @param filepath - path to file with domain describing
     * @return Domain structure
     * @throws std::runtime_error - when the file cannot been opened or the data in the file contains errors
     */
    auto Domain::readFromFile(const std::string& filepath)->Domain {
        return readDomainFromFile(filepath);
    }

    /**
     * @brief Make a dump of current domain includes additional comments
     * @return - string representation of dump file
     */
    auto Domain::dump() -> std::string {
        auto oss = ostringstream();

        oss << "************************\n";
        oss << "** Domain description **\n";
        oss << "************************\n\n";

        oss << format("// Count of coordinate lines for X and Y coordinates respectivly\n{} {}\n\n", this->Kx, this->Ky);

        oss << "// X-Y coordinate lines line-by-line\n";

        auto Xmax = findMaxSize(this->X) + 1;
        auto Ymax = findMaxSize(this->Y) + 1;
        for (size_t j = 0; j < this->Ky; j++) {
            for (size_t i = 0; i < this->Kx; i++) {
                size_t ind = i + this->Kx * j;
                oss << format("{0:>{2}} {1:<{3}}   ", this->X[ind], this->Y[ind], Xmax, Ymax);
            }
            oss << "\n";
        }
        oss << "\n\n";

        oss << format("// Count of subdomains\n{}\n\n", this->subdomains.size());

        oss << "// Subdomain describes in format `m[i] nxb[i] nxe[i] nyb[i] nye[i]` (*)\n";
        for (const auto& sd : this->subdomains) {
            oss << format("{:<3} {:>4} {:<4} {:>4} {:<4}\n", sd.materialNum, sd.xBeginNum, sd.xEndNum, sd.yBeginNum, sd.yEndNum);
        }
        oss << "\n";
        oss << "// (*)\n";
        oss << "// m[i]   - number of the material of the i'th subdomain\n";
        oss << "// nxb[i] - number X-coordinate line of begin of the i'th subdomain (starts from 1)\n";
        oss << "// nxe[i] - number X-coordinate line of end of the i'th subdomain (starts from 1)\n";
        oss << "// nyb[i] - number Y-coordinate line of begin of the i'th subdomain (starts from 1)\n";
        oss << "// nye[i] - number Y-coordinate line of end of the i'th subdomain (starts from 1)\n";
        oss << "\n\n";

        oss << "**********************\n";
        oss << "** Grid description **\n";
        oss << "**********************\n\n";

        oss << "// Number of subdivides and sparce coefficients for X coordinate lines in format like `nx[1] cx[1] ... nx[Kx-1] cx[Kx-1]`\n";
        for (size_t i = 0; i < this->nx.size(); i++) {
            oss << format("{} {}  ", this->nx.at(i), this->cx.at(i));
        }
        oss << "\n\n";
        oss << "// Number of subdivides and sparce coefficients for Y coordinate lines in format like `ny[1] cy[1] ... ny[Ky-1] cy[Ky-1]`\n";
        for (size_t i = 0; i < this->ny.size(); i++) {
            oss << format("{} {}  ", this->ny.at(i), this->cy.at(i));
        }
        oss << "\n\n";

        oss << format("// Number of additional subdivides for X and Y coordinate lines respectivly (0 is 'no additional subdivides') \n{} {}\n", this->splitX, this->splitY);

        return oss.str();
    }

    /**
     * @brief Make a dump of current domain without additional comments
     * @return - string representation of dump file
     */
    auto Domain::dumpNoComments() -> std::string {
        auto oss = ostringstream();

        oss << format("{} {} \n", this->Kx, this->Ky);

        for (size_t j = 0; j < this->Ky; j++) {
            for (size_t i = 0; i < this->Kx; i++) {
                size_t ind = i + this->Kx * j;
                oss << format("{} {} ", this->X[ind], this->Y[ind]);
            }
            oss << "\n";
        }
        oss << "\n";

        oss << format("{}\n", this->subdomains.size());

        for (const auto& sd : this->subdomains) {
            oss << format("{} {} {} {} {}\n", sd.materialNum, sd.xBeginNum, sd.xEndNum, sd.yBeginNum, sd.yEndNum);
        }
        oss << "\n";

        for (size_t i = 0; i < this->nx.size(); i++) {
            oss << format("{} {} ", this->nx.at(i), this->cx.at(i));
        }
        oss << "\n";
        for (size_t i = 0; i < this->ny.size(); i++) {
            oss << format("{} {} ", this->ny.at(i), this->cy.at(i));
        }
        oss << "\n\n";

        oss << format("{} {}\n", this->splitX, this->splitY);

        return oss.str();
    }

}

