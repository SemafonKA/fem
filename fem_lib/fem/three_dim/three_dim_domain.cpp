#include "three_dim_domain.h"

#include <sstream>
#include <string>
#include <format>
#include <algorithm>
#include <vector>
#include <fstream>

#include "../../logger.h"
#include "../../timer.h"

#undef max

using std::string;
using std::format;
using std::vector;
using std::ifstream;
using std::stringstream;
using std::istringstream;
using std::ostringstream;


[[nodiscard]]
static auto findMaxSize(const vector <double>& vec) -> size_t {
    size_t max = 0;
    for (const auto& el : vec) {
        max = std::max(max, format("{}", el).size());
    }
    return max;
}

[[nodiscard]]
static inline auto isComment(char c) -> bool {
    if (c == '*' || c == '/') {
        return true;
    }

    return false;
}


enum class States {
    coordinates_count,
    coordinates,
    coord_z_count,
    coord_z,
    subdomains_count,
    subdomains,
    edges_count,
    edges,
    subdivision_x,
    subdivision_y,
    splits,
    END
};

[[nodiscard]]
static auto stateName(States state) -> string {
    switch (state) {
        using enum States;
    case coordinates_count:
        return "coordinates_count";
    case coordinates:
        return "coordinates";

    case coord_z_count:
        return "coord_z_count";
    case coord_z:
        return "coord_z";

    case subdomains_count:
        return "subdomains_count";
    case subdomains:
        return "subdomains";

    case edges_count:
        return "edges_count";
    case edges:
        return "edges";

    case subdivision_x:
        return "subdivision_x";
    case subdivision_y:
        return "subdivision_y";

    case splits:
        return "splits";

    case END:
        return "END";

    default:
        throw std::runtime_error("Not implemented");
    }
}


class Reader {

public:
    auto readFromStr(const string& in) -> fem::three_dim::Domain {
        domain = fem::three_dim::Domain{};
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
                using enum States;
            case coordinates_count:
                success = readCoordinatesCount(token);
                break;

            case coordinates:
                success = readCoordinates(token);
                break;

            case coord_z_count:
                success = readCoordzCount(token);
                break;

            case coord_z:
                success = readCoordZ(token);
                break;

            case subdomains_count:
                success = readSubdomainCount(token);
                break;

            case subdomains:
                success = readSubdomains(token);
                break;

            case edges_count:
                success = readEdgesCount(token);
                break;

            case edges:
                success = readEdges(token);
                break;

            case subdivision_x:
                success = readSubdivX(token);
                break;

            case subdivision_y:
                success = readSubdivY(token);
                break;

            case splits:
                success = readSplits(token);
                break;

            case END:
                i = std::numeric_limits<size_t>::max();
                break;

            default:
                throw std::runtime_error("Not implemented branch");
            }

            if (!success) {
                string error = format("[line {}, char {}] - {}", lineNumber, charNumber, this->errorMsg);
                throw std::runtime_error(error);
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

    fem::three_dim::Domain domain = {};

    string errorMsg = ""; /// Contains error description if something goes wrong

    size_t edges_size = 0;

    /**
     * @brief Read next token from [in] start by [position]. Token is a "word" separated by whitespace, line break or EOF
     * @param in - string with content
     * @param position - start position to extract token
     * @return extracted token
     */
    string readToken(const string& in, size_t position) const noexcept {
        auto buf = ostringstream();
        size_t pos = position;
        while (pos < in.size() && in[pos] != ' ' && in[pos] != '\n') {
            buf << in[pos];
            pos++;
        }
        return buf.str();
    }

    /**
     * @brief Only set [errorMsg] and returns false
     * @param errText - text to set as error to
     * @return false as always
     */
    [[nodiscard]]
    inline bool error(const string& errText) noexcept {
        this->errorMsg = errText;
        return false;
    }

    bool readCoordinatesCount(const string& token) {
        size_t k = 0;

        if (auto buf = istringstream(token);  buf >> k) {
            if (k == 0) {
                if (domain.Kx != 0) {
                    return error("Ky value readed as zero, but must be greater than zero");
                }
                return error("Kx value readed as zero, but must be greater than zero");
            }

            if (domain.Kx != 0) {
                domain.Ky = k;
                state = States::coordinates;
            }
            else {
                domain.Kx = k;
            }
            return true;
        }

        if (domain.Kx != 0) {
            return error(format("Error while reading Ky number: an unsigned was expected, but '{}' was received", token));
        }
        return error(format("Error while reading Kx number: an unsigned was expected, but '{}' was received", token));
    }

    bool readCoordinates(const string& token) {
        double val = 0;

        if (auto buf = istringstream(token); buf >> val) {
            if (domain.X.size() > domain.Y.size()) {
                domain.Y.push_back(val);
            }
            else {
                domain.X.push_back(val);
            }
            if (domain.Y.size() == domain.Kx * domain.Ky) {
                state = States::coord_z_count;
            }
            return true;
        }

        return error(format("Error while reading coordinate lines: a float was expected, but '{}' was received", token));
    }

    bool readCoordzCount(const string& token) {
        size_t k = 0;

        if (auto buf = istringstream(token);  buf >> k) {
            if (k == 0) {
                return error("Kz value readed as zero, but must be greater than zero");
            }

            domain.Kz = k;
            state = States::coord_z;

            return true;
        }

        return error(format("Error while reading Kz number: an unsigned was expected, but '{}' was received", token));
    }

    bool readCoordZ(const string& token) {
        double val = 0;

        if (auto buf = istringstream(token); buf >> val) {
            domain.Z.push_back(val);

            if (domain.Z.size() == domain.Kz) {
                state = States::subdomains_count;
            }

            return true;
        }

        return error(format("Error while reading Z coordinates: a float was expected, but '{}' was received", token));
    }

    bool readSubdomainCount(const string& token) {
        size_t size = 0;

        if (auto buf = istringstream(token); buf >> size) {
            if (size == 0) {
                return error("Mx value readed as zero, but must be greater than zero");
            }

            domain.subdomains.resize(size);
            state = States::subdomains;
            return true;
        }

        return error(format("Error while reading Mx vakue: an unsigned was expected, but '{}' was received", token));
    }

    bool readSubdomains(const string& token) {
        static auto subdomains = vector<fem::three_dim::Subdomain>();
        static auto sd = fem::three_dim::Subdomain();
        static auto num = 0;

        size_t k = 0;
        if (auto buf = istringstream(token); buf >> k) {
            if (num == 0) {
                if (k == 0) {
                    return error("Error while reading material number: value must be greater than zero");
                }

                sd.materialNum = k;
                num++;
            }
            else {
                if (k == 0) {
                    return error("Error while reading subdomain coordinate line number: number must be greater than zero");
                }

                if (num == 1) {
                    if (k > domain.Kx) {
                        return error(format("Error while reading nxb number: number less than or equivalent to {} was expected, but {} was received", domain.Kx, k));
                    }

                    sd.xBeginNum = k;
                    num++;
                }
                else if (num == 2) {
                    if (k > domain.Kx) {
                        return error(format("Error while reading nxe number: number less than or equivalent to {} was expected, but {} was received", domain.Kx, k));
                    }
                    if (k <= sd.xBeginNum) {
                        return error(format("Error while reading nxe number: nxe number ({}) must be greater than nxb number ({})", k, sd.xBeginNum));
                    }

                    sd.xEndNum = k;
                    num++;
                }
                else if (num == 3) {
                    if (k > domain.Ky) {
                        return error(format("Error while reading nyb number: number less than or equivalent to {} was expected, but {} was received", domain.Ky, k));
                    }

                    sd.yBeginNum = k;
                    num++;
                }
                else if (num == 4) {
                    if (k > domain.Ky) {
                        return error(format("Error while reading nye number: number less than or equivalent to {} was expected, but {} was received", domain.Ky, k));
                    }
                    if (k <= sd.yBeginNum) {
                        return error(format("Error while reading nye number: nye number ({}) must be greater than nyb number ({})", k, sd.yBeginNum));
                    }

                    sd.yEndNum = k;
                    num++;
                }
                else if (num == 5) {
                    if (k > domain.Kz) {
                        return error(format("Error while reading nzb number: number less than or equivalent to {} was expected, but {} was received", domain.Kz, k));
                    }

                    sd.zBeginNum = k;
                    num++;
                }
                else {
                    if (k > domain.Kz) {
                        return error(format("Error while reading nze number: number less than or equivalent to {} was expected, but {} was received", domain.Kz, k));
                    }
                    if (k <= sd.zBeginNum) {
                        return error(format("Error while reading nze number: nze number ({}) must be greater than nzb number ({})", k, sd.zBeginNum));
                    }
                    sd.zEndNum = k;
                    num = 0;

                    subdomains.push_back(sd);
                    if (subdomains.size() == domain.subdomains.size()) {
                        domain.subdomains = std::move(subdomains);
                        state = States::edges_count;
                    }
                }
            }
            return true;
        }

        return error(format("Error while reading subdomains: unsigned was expected, but `{}` received", token));
    }

    bool readEdgesCount(const string& token) {
        size_t size = 0;

        if (auto buf = istringstream(token); buf >> size) {
            if (size == 0) {
                return error("Cx value readed as zero, but must be greater than zero");
            }

            edges_size = size;
            state = States::edges;
            return true;
        }

        return error(format("Error while reading Cx vakue: an unsigned was expected, but '{}' was received", token));
    }

    bool readEdges(const string& token) {
        static auto s1_edges = vector<fem::three_dim::S1Edge>();
        static auto s2_edges = vector<fem::three_dim::S2Edge>();
        static size_t num = 0;
        static size_t sType = 1;

        static size_t nxb = 0;
        static size_t nxe = 0;
        static size_t nyb = 0;
        static size_t nye = 0;
        static size_t nzb = 0;
        static size_t nze = 0;
        static size_t nc = 0;

        size_t k = 0;
        if (auto buf = istringstream(token); buf >> k) {
            if (num == 0) {
                if (k != 1 && k != 2) {
                    return error("Error while reading edge type: value must be '1' or '2'");
                }

                sType = k;
                num++;
            }
            else {
                if (k == 0) {
                    return error("Error while reading edge coordinate line number: number must be greater than zero");
                }

                if (num == 1) {
                    if (k > domain.Kx) {
                        return error(format("Error while reading nxb number: number less than or equivalent to {} was expected, but {} was received", domain.Kx, k));
                    }

                    nxb = k;
                    num++;
                }
                else if (num == 2) {
                    if (k > domain.Kx) {
                        return error(format("Error while reading nxe number: number less than or equivalent to {} was expected, but {} was received", domain.Kx, k));
                    }
                    if (k < nxb) {
                        return error(format("Error while reading nxe number: nxe number ({}) must be greater than nxb number ({})", k, nxb));
                    }

                    nxe = k;
                    num++;
                }
                else if (num == 3) {
                    if (k > domain.Ky) {
                        return error(format("Error while reading nyb number: number less than or equivalent to {} was expected, but {} was received", domain.Ky, k));
                    }

                    nyb = k;
                    num++;
                }
                else if (num == 4) {
                    if (k > domain.Ky) {
                        return error(format("Error while reading nye number: number less than or equivalent to {} was expected, but {} was received", domain.Ky, k));
                    }
                    if (k < nyb) {
                        return error(format("Error while reading nye number: nye number ({}) must be greater than nyb number ({})", k, nyb));
                    }

                    nye = k;
                    num++;
                }
                else if (num == 5) {
                    if (k > domain.Kz) {
                        return error(format("Error while reading nzb number: number less than or equivalent to {} was expected, but {} was received", domain.Kz, k));
                    }
                    nzb = k;
                    num++;
                }
                else if (num == 6) {
                    if (k > domain.Kz) {
                        return error(format("Error while reading nze number: number less than or equivalent to {} was expected, but {} was received", domain.Kz, k));
                    }
                    if (k < nzb) {
                        return error(format("Error while reading nze number: nze number ({}) must be greater than nzb number ({})", k, nzb));
                    }

                    nze = k;
                    num++;
                }
                else {
                    nc = k;
                    num = 0;

                    // Checks that edge condition setted up across one line and not many lines
                    if (nxb != nxe && nyb != nye && nzb != nze) {
                        return error(format(
                            "Error while reading edge condition: edge condition must be setted throw two coordinate lines.\n"
                            "                                        If [nxb] != [nxe] and [nyb] != [nye], then [nzb] must be equal to [nze]. \n"
                            "                                        If [nxb] != [nxe] and [nzb] != [nze], then [nyb] must be equal to [nye]. \n"
                            "                                        If [nyb] != [nye] and [nzb] != [nze], then [nxb] must be equal to [nxe]. \n"
                        ));
                    }

                    if (sType == 1) {
                        auto s1 = fem::three_dim::S1Edge();
                        s1.xBeginNum = nxb;
                        s1.xEndNum = nxe;

                        s1.yBeginNum = nyb;
                        s1.yEndNum = nye;

                        s1.zBeginNum = nzb;
                        s1.zEndNum = nze;
                        s1.funcNum = nc;

                        s1_edges.push_back(s1);
                    }
                    else {
                        auto s2 = fem::three_dim::S2Edge();
                        s2.xBeginNum = nxb;
                        s2.xEndNum = nxe;

                        s2.yBeginNum = nyb;
                        s2.yEndNum = nye;

                        s2.zBeginNum = nzb;
                        s2.zEndNum = nze;
                        s2.funcNum = nc;

                        s2_edges.push_back(s2);
                    }

                    if (s1_edges.size() + s2_edges.size() == edges_size) {
                        domain.s1_edges = std::move(s1_edges);
                        domain.s2_edges = std::move(s2_edges);

                        state = States::subdivision_x;
                    }
                }
            }
            return true;
        }

        return error(format("Error while reading subdomains: unsigned was expected, but `{}` received", token));
    }

    bool readSubdivX(const string& token) {
        auto buf = istringstream(token);
        if (domain.nx.size() > domain.cx.size()) {
            if (double val; buf >> val) {
                if (val <= 0.0) {
                    return error(format("Error while reading X subdivisions: sparse coef must be greater than 0, but {} was received", val));
                }

                domain.cx.push_back(val);
                if (domain.cx.size() == domain.Kx - 1) {
                    state = States::subdivision_y;
                }
                return true;
            }

            return error(format("Error while reading X subdivisions: float was expected, but `{}` received", token));
        }
        else {
            if (size_t k; buf >> k) {
                if (k == 0) {
                    return error("Error while reading X subdivisions: number of subdivisions must be greater than zero");
                }

                domain.nx.push_back(k);
                return true;
            }

            return error(format("Error while reading X subdivisions: unsigned was expected, but `{}` received", token));
        }
    }

    bool readSubdivY(const string& token) {
        auto buf = stringstream(token);
        if (domain.ny.size() > domain.cy.size()) {
            if (double val; buf >> val) {
                if (val <= 0.0) {
                    return error(format("Error while reading Y subdivisions: sparse coef must be greater than 0, but {} was received", val));
                }

                domain.cy.push_back(val);
                if (domain.cy.size() == domain.Ky - 1) {
                    state = States::splits;
                }
                return true;
            }

            return error(format("Error while reading Y subdivisions: float was expected, but `{}` received", token));
        }
        else {
            if (size_t k; buf >> k) {
                if (k == 0) {
                    return error("Error while reading Y subdivisions: number of subdivisions must be greater than zero");
                }

                domain.ny.push_back(k);
                return true;
            }
        }

        return error(format("Error while reading Y subdivisions: unsigned was expected, but `{}` received", token));
    }

    bool readSplits(const string& token) {
        static auto num = 0;

        size_t k = 0;
        if (auto buf = istringstream(token); buf >> k) {
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

        return error(format("Error while reading splits: unsigned was expected, but `{}` received", token));
    }
};

[[nodiscard]]
static auto readDomainFromFile(const std::string& filepath) -> fem::three_dim::Domain {
    auto file = std::ifstream(filepath);

    if (!file.is_open()) {
        const auto errorStr = format("in function `readDomainFromFile`: Error opening file \"{}\".", filepath);
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
        logger::debug(format("in function `readDomainFromFile`: {}", e.what()));
        throw std::runtime_error(e.what());
    }
}

namespace fem::three_dim {

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

        oss << format("// Count of coordinate lines by Z axis\n{}\n\n", this->Kz);
        auto Zmax = findMaxSize(this->Z) + 1;
        for (size_t i = 0; i < this->Z.size(); i++) {
            oss << format("{0:>{1}}   ", this->Z[i], Zmax);
        }
        oss << "\n\n";

        oss << format("// Count of subdomains\n{}\n\n", this->subdomains.size());

        oss << "// Subdomain describes in format `m[i] nxb[i] nxe[i] nyb[i] nye[i] nzb[i] nze[i]` (*)\n";
        for (const auto& sd : this->subdomains) {
            oss << format("{:<3} {:>4} {:<4} {:>4} {:<4} {:<4} {:<4}\n", sd.materialNum, sd.xBeginNum, sd.xEndNum, sd.yBeginNum, sd.yEndNum, sd.zBeginNum, sd.zEndNum);
        }
        oss << "\n";
        oss << "// (*)\n";
        oss << "// m[i]   - number of the material of the i'th subdomain\n";
        oss << "// nxb[i] - number X-coordinate line of begin of the i'th subdomain (starts from 1)\n";
        oss << "// nxe[i] - number X-coordinate line of end of the i'th subdomain (starts from 1)\n";
        oss << "// nyb[i] - number Y-coordinate line of begin of the i'th subdomain (starts from 1)\n";
        oss << "// nye[i] - number Y-coordinate line of end of the i'th subdomain (starts from 1)\n";
        oss << "// nzb[i] - number Z-coordinate line of begin of the i'th subdomain (starts from 1)\n";
        oss << "// nze[i] - number Z-coordinate line of end of the i'th subdomain (starts from 1)\n";
        oss << "\n\n";

        oss << format("// Count of edge conditions\n{}\n\n", this->s1_edges.size() + this->s2_edges.size());
        oss << "// Edge conditions describes in format `t[i] nxb[i] nxe[i] nyb[i] nye[i] nzb[i] nze[i] nc[i]` (**)\n";

        for (const auto& s1 : this->s1_edges) {
            oss << format("{:<3} {:>4} {:<4} {:>4} {:<4} {:<4} {:<4} {:>3}\n", 1, s1.xBeginNum, s1.xEndNum, s1.yBeginNum, s1.yEndNum, s1.zBeginNum, s1.zEndNum, s1.funcNum);
        }

        for (const auto& s2 : this->s2_edges) {
            oss << format("{:<3} {:>4} {:<4} {:>4} {:<4} {:>4} {:<4} {:>3}\n", 2, s2.xBeginNum, s2.xEndNum, s2.yBeginNum, s2.yEndNum, s2.zBeginNum, s2.zEndNum, s2.funcNum);
        }
        oss << "\n";
        oss << "// (*)\n";
        oss << "// t[i]   - type of the edge condition (1 or 2)\n";
        oss << "// nxb[i] - number X-coordinate line of begin of the i'th edge (starts from 1)\n";
        oss << "// nxe[i] - number X-coordinate line of end of the i'th edge (starts from 1)\n";
        oss << "// nyb[i] - number Y-coordinate line of begin of the i'th edge (starts from 1)\n";
        oss << "// nye[i] - number Y-coordinate line of end of the i'th edge (starts from 1)\n";
        oss << "// nzb[i] - number Z-coordinate line of begin of the i'th edge (starts from 1)\n";
        oss << "// nze[i] - number Z-coordinate line of end of the i'th edge (starts from 1)\n";
        oss << "// nc[i]  - number of the edge condition function (related to function nums in program)\n";
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

        oss << this->Kz << "\n";
        for (const auto& v : this->Z) {
            oss << v << " ";
        }
        oss << "\n\n";

        oss << format("{}\n", this->subdomains.size());

        for (const auto& sd : this->subdomains) {
            oss << format("{} {} {} {} {} {} {}\n", sd.materialNum, sd.xBeginNum, sd.xEndNum, sd.yBeginNum, sd.yEndNum, sd.zBeginNum, sd.zEndNum);
        }
        oss << "\n";

        for (const auto& s1 : this->s1_edges) {
            oss << format("{} {} {} {} {} {} {} {}\n", 1, s1.xBeginNum, s1.xEndNum, s1.yBeginNum, s1.yEndNum, s1.zBeginNum, s1.zEndNum, s1.funcNum);
        }
        for (const auto& s2 : this->s2_edges) {
            oss << format("{} {} {} {} {} {} {} {}\n", 2, s2.xBeginNum, s2.xEndNum, s2.yBeginNum, s2.yEndNum, s2.zBeginNum, s2.zEndNum, s2.funcNum);
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
