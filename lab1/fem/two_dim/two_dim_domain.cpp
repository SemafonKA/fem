#include "two_dim_domain.h"

#include "two_dim_domain_reader.h"

namespace fem::two_dim {

    /**
     * @brief Read domain structure like described in example from file (ignoring comments)
     * @param filepath - path to file with domain describing
     * @return Domain structure
     * @throws std::runtime_error - when the file cannot been opened or the data in the file contains errors
     */
    auto fem::two_dim::Domain::readFromFile(const std::string& filepath)->fem::two_dim::Domain {
        return readDomainFromFile(filepath);
    }

}
