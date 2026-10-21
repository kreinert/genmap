#include <gtest/gtest.h>

#include <algorithm>
#include <chrono>

#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>

static constexpr bool outputProgress = false;

#include "../src/common.hpp"
#include "../src/algo.hpp"

using namespace seqan;

std::mt19937_64 rng;

template <typename TSpec, typename TLengthSum, unsigned LEVELS, unsigned WORDS_PER_BLOCK>
unsigned GemMapFastFMIndexConfig<TSpec, TLengthSum, LEVELS, WORDS_PER_BLOCK>::SAMPLING = 10;

template <typename TChar, typename TSpec, typename TRng>
void randomText(String<TChar, TSpec> & string, TRng & rng, uint64_t const length)
{
    std::uniform_int_distribution<typename TRng::result_type> distr(0, ValueSize<TChar>::VALUE - 1);
    resize(string, length);
    for (uint64_t i = 0; i < length; ++i)
        string[i] = TChar(distr(rng));
}

template <typename TIndexIt, typename TNeedle, typename TNeedleIt, typename TThreshold, typename TDistance>
inline void
_trivialBacktracking(TIndexIt indexIt, TNeedle const & needle, TNeedleIt needleIt,
                     TThreshold errors, TThreshold threshold, uint64_t & frequency, TDistance)
{
    constexpr bool isDna5 = std::is_same<typename Value<TNeedle>::Type, Dna5>::value;

    // Exact case.
    if (errors == threshold)
    {
        while (!atEnd(needleIt, needle))
        {
            if ((isDna5 && (*needleIt == Dna5('N'))) || !goDown(indexIt, *needleIt))
                break;
            ++needleIt;
        }
        if (atEnd(needleIt, needle))
            frequency += countOccurrences(indexIt);
    }
    // Approximate case.
    else if (errors < threshold)
    {
        // Base case.
        if (atEnd(needleIt, needle))
        {
            frequency += countOccurrences(indexIt);
        }
        // Recursive case.
        else
        {
            // Insertion.
            SEQAN_IF_CONSTEXPR (IsSameType<TDistance, EditDistance>::VALUE)
            {
                _trivialBacktracking(indexIt, needle, needleIt + 1,
                                     static_cast<TThreshold>(errors + 1), threshold, frequency, TDistance());
            }

            if (goDown(indexIt))
            {
                do
                {
                    // Mismatch.
                    TThreshold delta = (isDna5 && *needleIt == Dna5('N')) || !ordEqual(parentEdgeLabel(indexIt), *needleIt);
                    _trivialBacktracking(indexIt, needle, needleIt + 1,
                                         static_cast<TThreshold>(errors + delta), threshold, frequency, TDistance());

                    // Deletion.
                    SEQAN_IF_CONSTEXPR (IsSameType<TDistance, EditDistance>::VALUE)
                    {
                        _trivialBacktracking(indexIt, needle, needleIt,
                                             static_cast<TThreshold>(errors + 1), threshold, frequency, TDistance());
                    }
                }
                while (goRight(indexIt));
            }
        }
    }
}

template <typename TIndex, typename TNeedle, typename TThreshold, typename TDistance>
uint64_t trivialBacktracking(TIndex & index, TNeedle const & needle, TThreshold threshold, TDistance const & /*tag*/)
{
    typedef typename Iterator<TIndex, TopDown<> >::Type       TIndexIt;
    typedef typename Iterator<TNeedle const, Standard>::Type  TNeedleIt;

    TIndexIt indexIt(index);
    TNeedleIt needleIt = begin(needle, Standard());
    TThreshold errors = 0;

    uint64_t frequency = 0;
    _trivialBacktracking(indexIt.revIter, needle, needleIt, errors, threshold, frequency, TDistance());
    return frequency;
}

template <typename TDistance, typename TChar, typename TIndex, typename TContainer>
inline void computeMappabilityTrivial(TIndex & index, TContainer & c, SearchParams const & searchParams, unsigned const errors)
{
    using value_type = typename TContainer::value_type;

    constexpr uint64_t max_val = std::numeric_limits<typename TContainer::value_type>::max();

    auto const & text = indexText(index);

    uint64_t global_pos = 0;
    for (uint64_t seq = 0; seq < length(text); ++seq)
    {
        for (uint64_t i = 0; i < length(text[seq]) - searchParams.length + 1; ++i, ++global_pos)
        {
            auto const & needle = infix(text[seq], i, i + searchParams.length);

            value_type hits = std::min(trivialBacktracking(index, needle, errors, TDistance()), max_val);
            if (searchParams.revCompl && hits < max_val)
            {
                String<TChar> needleRevCompl(needle);
                reverseComplement(needleRevCompl);
                hits = std::min(hits + trivialBacktracking(index, needleRevCompl, errors, TDistance()), max_val);
            }
            c[global_pos] = hits;
        }
        global_pos += searchParams.length - 1;
    }
}

template <typename TChar, typename TDistance, unsigned errors>
void test(uint64_t const nbrChromosomes, uint64_t const lengthChromosomes, uint64_t const iterations)
{
    using TIndexConfig = TBiIndexConfig<TGemMapFastFMIndexConfig<uint32_t>>;

    for (uint64_t it = 0; it < iterations; ++it)
    {
        typedef StringSet<String<TChar>, Owner<ConcatDirect<> > > TGenome;
        TGenome genome;

        // TODO: replace with stringSetLimits
        StringSet<uint64_t> chromLengths, chromCumLengths; // needed for localization and reset

        uint64_t cumLength = 0;
        appendValue(chromCumLengths, 0);
        for (uint64_t ss = 0; ss < nbrChromosomes; ++ss)
        {
            String<TChar> chr;
            randomText(chr, rng, lengthChromosomes);
            appendValue(genome, chr);
            appendValue(chromLengths, lengthChromosomes);
            cumLength += lengthChromosomes;
            appendValue(chromCumLengths, cumLength);
        }
        // auto const chromLengths = stringSetLimits(genome);

        Index<TGenome, TIndexConfig> index(genome);
        indexCreate(index, FibreSALF());
        auto const & text = indexText(index).concat;

        uint64_t const totalLength = seqan::length(text);
        std::vector<uint8_t> frequencyGenMap(totalLength), frequencyTrivial(totalLength);

        uint64_t const minK = errors + 1 + (errors >= 2);

        for (uint64_t k = minK; k <= 8; ++k)
        {
            SearchParams searchParams;
            searchParams.length = k;
            searchParams.threads = omp_get_num_threads();
            searchParams.revCompl = rng() % 2;
            searchParams.excludePseudo = false;

            frequencyTrivial.assign(totalLength, 0);
            computeMappabilityTrivial<TDistance, TChar>(index, frequencyTrivial, searchParams, errors);

            // iterate over all possible overlap values
            for (uint64_t overlap = minK; overlap <= k; ++overlap)
            {
                searchParams.overlap = overlap;
                // std::cout << "E: " << errors << ", K: " << k << ", O: " << overlap << std::endl;

                using TLocation = Pair<uint16_t, uint32_t>;
                std::map<TLocation, std::pair<std::vector<TLocation>, std::vector<TLocation> > > locations;
                std::vector<uint16_t> mappingSeqIdFile(0);

                frequencyGenMap.assign(totalLength, 0);
                std::vector<std::pair<uint64_t, uint64_t> > intervals;
                bool completeSameKmers;
                computeMappability<errors>(index, text, frequencyGenMap, searchParams, false /*dir*/, chromLengths, chromCumLengths,
                                           locations, mappingSeqIdFile, intervals, completeSameKmers, 1/*currentFileNo*/, 1/*totalFileNo*/, false /*csvComputation*/);

                EXPECT_EQ(frequencyTrivial, frequencyGenMap);
                // if (frequencyTrivial != frequencyGenMap)
                // {
                //     std::cerr << "K: " << k << ", Overlap: " << overlap << '\n';
                //     for (uint64_t ss = 0; ss < nbrChromosomes; ++ss)
                //         std::cerr << genome[ss] << '\n';
                //     std::copy(frequencyTrivial.begin(), frequencyTrivial.end(), std::ostream_iterator<int>(std::cerr, " "));
                //     std::cerr << '\n';
                //     std::copy(frequencyGenMap.begin(), frequencyGenMap.end(), std::ostream_iterator<int>(std::cerr, " "));
                //     std::cerr << '\n';
                //     exit(1);
                // }
            }
        }
    }
}

template <typename TIndex, typename TNeedle, typename TLimits>
void trivialEditCollect(Iter<TIndex, VSTree<TopDown<> > > it, TNeedle const & needle, uint32_t const needlePos,
                        unsigned const errors, unsigned const threshold, std::vector<uint64_t> & starts, TLimits const & limits)
{
    constexpr bool isDna5 = std::is_same<typename Value<TNeedle>::Type, Dna5>::value;
    if (errors > threshold)
        return;

    if (needlePos == length(needle))
    {
        for (auto const & occ : getOccurrences(it))
            starts.push_back(posGlobalize(occ, limits));
        return;
    }

    // Insertion in the pattern: the text substring is shorter.
    if (errors < threshold)
        trivialEditCollect(it, needle, needlePos + 1, errors + 1, threshold, starts, limits);

    if (goDown(it, Rev()))
    {
        do
        {
            unsigned const delta = (!ordEqual(parentEdgeLabel(it, Rev()), needle[needlePos])
                                    || (isDna5 && needle[needlePos] == Dna5('N'))) ? 1u : 0u;
            if (errors + delta <= threshold)
                trivialEditCollect(it, needle, needlePos + 1, errors + delta, threshold, starts, limits);
            // Deletion in the pattern: the text substring is longer.
            if (errors + 1 <= threshold)
                trivialEditCollect(it, needle, needlePos, errors + 1, threshold, starts, limits);
        } while (goRight(it, Rev()));
    }
}

template <typename TIndex, typename TNeedle>
uint64_t trivialEditFrequency(TIndex & index, TNeedle const & needle, unsigned const errors)
{
    std::vector<uint64_t> starts;
    trivialEditCollect(Iter<TIndex, VSTree<TopDown<> > >(index), needle, 0, 0, errors, starts, stringSetLimits(indexText(index)));
    std::sort(starts.begin(), starts.end());
    starts.erase(std::unique(starts.begin(), starts.end()), starts.end());
    return starts.size();
}

template <typename TChar, unsigned errors>
void testEdit(uint64_t const nbrChromosomes, uint64_t const lengthChromosomes, uint64_t const kFrom, uint64_t const kTo, bool const revCompl)
{
    using TIndexConfig = TBiIndexConfig<TGemMapFastFMIndexConfig<uint32_t>>;
    typedef StringSet<String<TChar>, Owner<ConcatDirect<> > > TGenome;
    TGenome genome;

    StringSet<uint64_t> chromLengths, chromCumLengths;
    uint64_t cumLength = 0;
    appendValue(chromCumLengths, 0);
    for (uint64_t ss = 0; ss < nbrChromosomes; ++ss)
    {
        String<TChar> chr;
        randomText(chr, rng, lengthChromosomes);
        appendValue(genome, chr);
        appendValue(chromLengths, lengthChromosomes);
        cumLength += lengthChromosomes;
        appendValue(chromCumLengths, cumLength);
    }

    Index<TGenome, TIndexConfig> index(genome);
    indexCreate(index, FibreSALF());
    auto const & text = indexText(index).concat;
    uint64_t const totalLength = seqan::length(text);

    for (uint64_t k = kFrom; k <= kTo; ++k)
    {
        SearchParams searchParams;
        searchParams.length = k;
        searchParams.overlap = 0;
        searchParams.threads = 2;
        searchParams.revCompl = revCompl;
        searchParams.excludePseudo = false;

        std::vector<uint8_t> frequencyEdit(totalLength, 0), frequencyOracle(totalLength, 0);
        using TLocation = Pair<uint16_t, uint32_t>;
        std::map<TLocation, std::pair<std::vector<TLocation>, std::vector<TLocation> > > locations;
        std::vector<uint16_t> mappingSeqIdFile(0);
        std::vector<std::pair<uint64_t, uint64_t> > intervals;
        bool completeSameKmers = false;
        computeMappabilityEdit<errors>(index, text, frequencyEdit, searchParams, false, chromLengths, chromCumLengths,
                                       locations, mappingSeqIdFile, intervals, completeSameKmers, 1, 1, false);

        constexpr uint64_t max_val = std::numeric_limits<uint8_t>::max();
        uint64_t globalPos = 0;
        for (uint64_t seq = 0; seq < length(genome); ++seq)
        {
            for (uint64_t i = 0; i + k <= length(genome[seq]); ++i, ++globalPos)
            {
                auto const needle = infix(genome[seq], i, i + k);
                uint64_t hits = trivialEditFrequency(index, needle, errors);
                if (revCompl && hits < max_val)
                {
                    String<TChar> needleRevCompl(needle);
                    reverseComplement(needleRevCompl);
                    hits = std::min(hits + trivialEditFrequency(index, needleRevCompl, errors), max_val);
                }
                frequencyOracle[globalPos] = static_cast<uint8_t>(std::min(hits, max_val));
            }
            globalPos += k - 1;
        }

        EXPECT_EQ(frequencyOracle, frequencyEdit);
    }
}

TEST(GenMapAlgo, edit_0_matches_hamming)
{
    using TIndexConfig = TBiIndexConfig<TGemMapFastFMIndexConfig<uint32_t>>;
    typedef StringSet<String<Dna>, Owner<ConcatDirect<> > > TGenome;
    TGenome genome;
    StringSet<uint64_t> chromLengths, chromCumLengths;
    appendValue(chromCumLengths, 0);
    uint64_t cumLength = 0;
    for (uint64_t ss = 0; ss < 2; ++ss)
    {
        String<Dna> chr;
        randomText(chr, rng, 40);
        appendValue(genome, chr);
        appendValue(chromLengths, 40);
        cumLength += 40;
        appendValue(chromCumLengths, cumLength);
    }
    Index<TGenome, TIndexConfig> index(genome);
    indexCreate(index, FibreSALF());
    auto const & text = indexText(index).concat;

    SearchParams searchParams;
    searchParams.length = 5;
    searchParams.overlap = 5;
    searchParams.threads = 2;
    searchParams.revCompl = true;
    searchParams.excludePseudo = false;

    std::vector<uint8_t> hamming(length(text), 0), edit(length(text), 0);
    using TLocation = Pair<uint16_t, uint32_t>;
    std::map<TLocation, std::pair<std::vector<TLocation>, std::vector<TLocation> > > locations;
    std::vector<uint16_t> mappingSeqIdFile(0);
    std::vector<std::pair<uint64_t, uint64_t> > intervals;
    bool completeSameKmers = false;
    computeMappability<0>(index, text, hamming, searchParams, false, chromLengths, chromCumLengths,
                          locations, mappingSeqIdFile, intervals, completeSameKmers, 1, 1, false);
    computeMappabilityEdit<0>(index, text, edit, searchParams, false, chromLengths, chromCumLengths,
                              locations, mappingSeqIdFile, intervals, completeSameKmers, 1, 1, false);
    EXPECT_EQ(hamming, edit);
}

TEST(GenMapAlgo, edit_1_dna4)
{
    testEdit<Dna, 1>(2, 18, 2, 6, false);
    testEdit<Dna, 1>(2, 18, 2, 5, true);
}

TEST(GenMapAlgo, edit_2_dna4)
{
    testEdit<Dna, 2>(1, 14, 4, 6, false);
}

TEST(GenMapAlgo, edit_1_dna5)
{
    testEdit<Dna5, 1>(1, 12, 2, 5, true);
}

TEST(GenMapAlgo, exact_dna4)
{
    test<Dna, HammingDistance, 0>(3, 1000, 1);
}

TEST(GenMapAlgo, hamming_1_dna4)
{
    test<Dna, HammingDistance, 1>(3, 1000, 1);
}

TEST(GenMapAlgo, hamming_2_dna4)
{
    test<Dna, HammingDistance, 2>(3, 1000, 1);
}

TEST(GenMapAlgo, hamming_3_dna4)
{
    test<Dna, HammingDistance, 3>(3, 1000, 1);
}

TEST(GenMapAlgo, hamming_4_dna4)
{
    test<Dna, HammingDistance, 4>(3, 1000, 1);
}

TEST(GenMapAlgo, exact_dna5)
{
    test<Dna5, HammingDistance, 0>(3, 1000, 1);
}

TEST(GenMapAlgo, hamming_1_dna5)
{
    test<Dna5, HammingDistance, 1>(3, 1000, 1);
}

TEST(GenMapAlgo, hamming_2_dna5)
{
    test<Dna5, HammingDistance, 2>(3, 1000, 1);
}

TEST(GenMapAlgo, hamming_3_dna5)
{
    test<Dna5, HammingDistance, 3>(3, 1000, 1);
}

TEST(GenMapAlgo, hamming_4_dna5)
{
    test<Dna5, HammingDistance, 4>(3, 1000, 1);
}

int main(int argc, char ** argv)
{
    auto now = std::chrono::system_clock::now();
    auto seed = std::chrono::duration_cast<std::chrono::nanoseconds>(now.time_since_epoch()).count();
    rng.seed(seed);
    std::cout << "Seed: " << seed << '\n';

    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
