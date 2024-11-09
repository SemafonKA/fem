#include "domain_reader.h"

#include <string>
#include <sstream>
#include <fstream>
#include <format>

#include "../../logger.h"

#undef max

using std::format;
using std::vector;
using std::string;
using std::ifstream;
using std::istringstream;
using std::stringstream;


static inline
auto isComment(char c) -> bool {
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

namespace fem::two_dim {

    auto readDomainFromFile(const std::string& filepath) -> fem::two_dim::Domain {
        auto file = std::ifstream(filepath);

        if (!file.is_open()) {
            const auto errorStr = format("Error opening file \"{}\".", filepath);
            logger::error(errorStr);
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
            logger::error(e.what());
            throw e;
        }
    }

}
