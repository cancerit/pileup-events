#include <catch2/catch_test_macros.hpp>

#include "const.hpp"
#include "pileup.hpp"

TEST_CASE ("score single") {
    std::vector<int> res (24,
                          0); // init result arr, single position long
    count_params cp; // defaults fine
    AEVSettings aevst; // defaults fine (discard_overlaps=false)
    AlleleEventCounter aev (cp, res, aevst);

    BaseInfo b{HTS_NT_A, 0, 30}; // a very normal a

    aev._score_single (b, 0);
    CAPTURE (res);
    REQUIRE (res[FIELD_A] == 1);
    REQUIRE (res[FIELD_MAPQ] == 30);
    REQUIRE (res[FIELD_NOBS] == 1);
    int sum = 0;
    for (int x : res)
        sum += x;
    REQUIRE (sum == 32);  // rest 0
}
