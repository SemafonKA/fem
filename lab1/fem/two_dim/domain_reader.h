#pragma once
#include <string>

#include "domain.h"


namespace fem::two_dim {

    /**
     * @brief Read info about domain from file `filepath`
     * @param filepath - path to file
     * @return - readed domain
     */
    auto readDomainFromFile(const std::string& filepath) -> fem::two_dim::Domain;

}
