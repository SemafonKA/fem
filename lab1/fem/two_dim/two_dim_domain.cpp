#include "two_dim_domain.h"

#include <sstream>
#include <string>
#include <format>
#include <algorithm>
#include <vector>
#include <fstream>

#include "../../logger.h"
#include "../../timer.h"

#undef max

using std::ostringstream;
using std::string;
using std::format;
using std::vector;
using std::ifstream;
using std::istringstream;
using std::stringstream;


static auto findMaxSize(const vector <double>& vec) -> size_t {
    size_t max = 0;
    for (const auto& el : vec) {
        max = std::max(max, format("{}", el).size());
    }
    return max;
}

static inline auto isComment(char c) -> bool {
    if (c == '*' || c == '/') {
        return true;
    }

    return false;
}


enum class States {
    coordinates_count,
    coordinates,
    subdomains_count,
    subdomains,
    subdivision_x,
    subdivision_y,
    splits,
    END
};

string stateName(States state) {
    switch (state) {
    case States::coordinates_count:
        return "coordinates_count";

    case States::coordinates:
        return "coordinates";

    case States::subdomains_count:
        return "subdomains_count";

    case States::subdomains:
        return "subdomains";

    case States::subdivision_x:
        return "subdivision_x";

    case States::subdivision_y:
        return "subdivision_y";

    case States::splits:
        return "splits";

    case States::END:
        return "END";

    default:
        throw std::runtime_error("Not implemented");
    }
}


class Reader {

public:
    auto readFromStr(const string& in) -> fem::two_dim::Domain {
        domain = fem::two_dim::Domain{};
        state = States::coordinates_count;

        size_t lineNumber = 1;
        size_t charNumber = 1;
        bool comment = false;

        for (size_t i = 0; i < in.size(); i++) {
            if (in[i] == '\n') {
                comment = false;
                charNumber = 1;
                lineNumber++;
                continue;
            }

            if (in[i] == ' ') {
                charNumber++;
                continue;
            }

            if (comment || isComment(in[i])) {
                comment = true;
                continue;
            }

            auto token = readToken(in, i);
            bool success = false;

            switch (state) {
            case States::coordinates_count:
                success = readCoordinatesCount(token);
                if (!success) {
                    string error = format("[line {}, char {}] - Error when reading Kx Ky numbers", lineNumber, charNumber);
                    throw std::runtime_error(error);
                }
                break;

            case States::coordinates:
                success = readCoordinates(token);
                if (!success) {
                    string error = format("[line {}, char {}] - Error when reading coordinate lines", lineNumber, charNumber);
                    throw std::runtime_error(error);
                }
                break;

            case States::subdomains_count:
                success = readSubdomainCount(token);
                if (!success) {
                    string error = format("[line {}, char {}] - Error when reading No number", lineNumber, charNumber);
                    throw std::runtime_error(error);
                }
                break;

            case States::subdomains:
                success = readSubdomains(token);
                if (!success) {
                    string error = format("[line {}, char {}] - Error when reading subdomains", lineNumber, charNumber);
                    throw std::runtime_error(error);
                }
                break;

            case States::subdivision_x:
                success = readSubdivX(token);
                if (!success) {
                    string error = format("[line {}, char {}] - Error when reading X subdivisions", lineNumber, charNumber);
                    throw std::runtime_error(error);
                }
                break;

            case States::subdivision_y:
                success = readSubdivY(token);
                if (!success) {
                    string error = format("[line {}, char {}] - Error when reading Y subdivisions", lineNumber, charNumber);
                    throw std::runtime_error(error);
                }
                break;

            case States::splits:
                success = readSplits(token);
                if (!success) {
                    string error = format("[line {}, char {}] - Error when reading splits", lineNumber, charNumber);
                    throw std::runtime_error(error);
                }
                break;

            case States::END:
                i = std::numeric_limits<size_t>::max();
                break;

            default:
                throw std::runtime_error("Not implemented branch");
            }

            if (success) {
                charNumber += token.length();
                i += token.length() - 1;
            }
        }

        if (state != States::END) {
            string error = format("Not enought data. Last readins state is: {}", stateName(state));
            throw std::runtime_error(error);
        }

        return std::move(this->domain);
    }

private:
    States state = States::coordinates_count;

    fem::two_dim::Domain domain = {};

    string readToken(const string& in, size_t position) const noexcept {
        auto buf = stringstream();
        size_t pos = position;
        while (pos < in.size() && in[pos] != ' ' && in[pos] != '\n') {
            buf << in[pos];
            pos++;
        }
        return buf.str();
    }

    bool readCoordinatesCount(const string& token) {
        size_t k = 0;
        if (auto buf = stringstream(token);  buf >> k) {
            if (domain.Kx != 0) {
                domain.Ky = k;
                state = States::coordinates;
            }
            else {
                domain.Kx = k;
            }
            return true;
        }
        return false;
    }

    bool readCoordinates(const string& token) {
        double val = 0;

        if (auto buf = stringstream(token); buf >> val) {
            if (domain.X.size() > domain.Y.size()) {
                domain.Y.push_back(val);
            }
            else {
                domain.X.push_back(val);
            }
            if (domain.Y.size() == domain.Kx * domain.Ky) {
                state = States::subdomains_count;
            }
            return true;
        }
        return false;
    }

    bool readSubdomainCount(const string& token) {
        size_t size = 0;

        if (auto buf = stringstream(token); buf >> size) {
            domain.subdomains.resize(size);
            state = States::subdomains;
            return true;
        }
        return false;
    }

    bool readSubdomains(const string& token) {
        static auto subdomains = vector<fem::two_dim::Subdomain>();
        static auto sd = fem::two_dim::Subdomain();
        static auto num = 0;

        size_t k = 0;
        if (auto buf = stringstream(token); buf >> k) {
            if (num == 0) {
                sd.materialNum = k;
                num++;
            }
            else if (num == 1) {
                sd.xBeginNum = k;
                num++;
            }
            else if (num == 2) {
                sd.xEndNum = k;
                num++;
            }
            else if (num == 3) {
                sd.yBeginNum = k;
                num++;
            }
            else {
                sd.yEndNum = k;
                num = 0;

                subdomains.push_back(sd);
                if (subdomains.size() == domain.subdomains.size()) {
                    domain.subdomains = std::move(subdomains);
                    state = States::subdivision_x;
                }
            }
            return true;
        }

        return false;
    }

    bool readSubdivX(const string& token) {
        auto buf = stringstream(token);
        if (domain.nx.size() > domain.cx.size()) {
            if (double val; buf >> val) {
                domain.cx.push_back(val);
                if (domain.cx.size() == domain.Kx - 1) {
                    state = States::subdivision_y;
                }
                return true;
            }
        }
        else {
            if (size_t k; buf >> k) {
                domain.nx.push_back(k);
                return true;
            }
        }

        return false;
    }

    bool readSubdivY(const string& token) {
        auto buf = stringstream(token);
        if (domain.ny.size() > domain.cy.size()) {
            if (double val; buf >> val) {
                domain.cy.push_back(val);
                if (domain.cy.size() == domain.Ky - 1) {
                    state = States::splits;
                }
                return true;
            }
        }
        else {
            if (size_t k; buf >> k) {
                domain.ny.push_back(k);
                return true;
            }
        }

        return false;
    }

    bool readSplits(const string& token) {
        static auto num = 0;

        size_t k = 0;
        if (auto buf = stringstream(token); buf >> k) {
            if (num == 0) {
                domain.splitX = k;
                num = 1;
            }
            else {
                domain.splitY = k;
                num = 0;
                state = States::END;
            }
            return true;
        }
        return false;
    }
};

static auto readDomainFromFile(const std::string& filepath) -> fem::two_dim::Domain {
    auto file = std::ifstream(filepath);

    if (!file.is_open()) {
        const auto errorStr = format("Error opening file \"{}\".", filepath);
        logger::debug(errorStr);
        throw std::runtime_error(errorStr);
    }
    // Read entire file at once
    auto fileContent = string((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    file.close();

    auto reader = Reader();
    try {
        auto result = reader.readFromStr(fileContent);
        return result;
    }
    catch (std::exception& e) {
        logger::debug(e.what());
        throw std::runtime_error(e.what());
    }
}

namespace fem::two_dim {

    /**
     * @brief Read domain structure like described in example from file (ignoring comments)
     * @param filepath - path to file with domain describing
     * @return Domain structure
     * @throws std::runtime_error - when the file cannot been opened or the data in the file contains errors
     */
    auto Domain::readFromFile(const std::string& filepath)->Domain {
        auto timer = Timer();
        auto domain = readDomainFromFile(filepath);
        auto elapsed = timer.elapsedMilliseconds();
        logger::debug(format("Reading domain from file completed by {} ms", elapsed));

        return domain;
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
                oss << format("{:20.14e} {:20.14e} ", this->X[ind], this->Y[ind]);
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
            oss << format("{} {:20.14e} ", this->nx.at(i), this->cx.at(i));
        }
        oss << "\n";
        for (size_t i = 0; i < this->ny.size(); i++) {
            oss << format("{} {:20.14e} ", this->ny.at(i), this->cy.at(i));
        }
        oss << "\n\n";

        oss << format("{} {}\n", this->splitX, this->splitY);

        return oss.str();
    }

}

