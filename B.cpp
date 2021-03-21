#include <iostream>
#include <vector>
#include <string>

using std::cin;
using std::cout;
using std::vector;
using std::string;
using std::max;

class SuffixArrayWithLCP {
public:
    SuffixArrayWithLCP() = default;

    explicit SuffixArrayWithLCP(string input);

    size_t operator[](size_t pos) const;

    size_t GetLCP(size_t pos) const;

private:
    struct RankedSuffix {
        size_t index;
        int rank;
        int rank_next;
    };

    void BuildSuffixArray(vector<RankedSuffix> &suffixes);

    void BuildSuffix();

    void BuildLCP();

    static const size_t alphabet_size = 255;
    string string_;
    vector<size_t> suffix_array_;
    vector<size_t> lcp_array_;
};


const size_t SuffixArrayWithLCP::alphabet_size;

SuffixArrayWithLCP::SuffixArrayWithLCP(string input) : string_(move(input)) {
    BuildSuffix();
    BuildLCP();
}

size_t SuffixArrayWithLCP::operator[](size_t pos) const {
    return suffix_array_[pos];
}

size_t SuffixArrayWithLCP::GetLCP(size_t pos) const {
    return lcp_array_[pos];
}

void SuffixArrayWithLCP::BuildSuffix() {
    size_t size = string_.length();
    vector<RankedSuffix> suffixes(size);

    for (size_t i = 0; i < size; ++i) {
        suffixes[i].index = i;
    }

    for (size_t cycle_size = 2; cycle_size < 2 * size; cycle_size <<= 1) {
        if (cycle_size == 2) {
            for (size_t i = 0; i < size; ++i) {
                suffixes[i].rank = string_[i];
                suffixes[i].rank_next = i + 1 < size ? string_[i + 1] : -1;
            }
        } else {
            int prev_rank = 0;
            vector<size_t> index_to_suffix(size);
            for (size_t i = 0; i < size; ++i) {
                if (i == 0) {
                    prev_rank = suffixes[i].rank;
                    suffixes[i].rank = 0;
                } else {
                    if (suffixes[i].rank == prev_rank &&
                        suffixes[i].rank_next == suffixes[i - 1].rank_next) {
                        suffixes[i].rank = suffixes[i - 1].rank;
                    } else {
                        prev_rank = suffixes[i].rank;
                        suffixes[i].rank = suffixes[i - 1].rank + 1;
                    }
                    index_to_suffix[suffixes[i].index] = i;
                }
            }

            for (size_t i = 0; i < size; ++i) {
                size_t next_index = suffixes[i].index + cycle_size / 2;
                suffixes[i].rank_next = next_index < size
                                        ? suffixes[index_to_suffix[next_index]].rank
                                        : -1;
            }
        }

        BuildSuffixArray(suffixes);
    }
}

void SuffixArrayWithLCP::BuildSuffixArray(vector<RankedSuffix> &suffixes) {
    size_t size = string_.length();
    for (int rank_type = 0; rank_type < 2; ++rank_type) {
        vector<RankedSuffix> buffer(size);
        vector<size_t> counter(max<size_t>(size + 1, alphabet_size));
        for (int i = 0; i < size; ++i) {
            int value = rank_type ? suffixes[i].rank : suffixes[i].rank_next;
            ++counter[value + 1];
        }
        for (int i = 1; i < counter.size(); ++i) {
            counter[i] += counter[i - 1];
        }
        for (size_t i = size; i > 0; --i) {
            int value = rank_type ? suffixes[i - 1].rank : suffixes[i - 1].rank_next;
            buffer[counter[value + 1] - 1] = suffixes[i - 1];
            counter[value + 1]--;
        }
        for (size_t i = 0; i < size; i++) {
            suffixes[i] = buffer[i];
        }
    }

    suffix_array_.resize(size);
    for (int i = 0; i < size; ++i) {
        suffix_array_[i] = suffixes[i].index;
    }
}


void SuffixArrayWithLCP::BuildLCP() {
    const size_t size = suffix_array_.size();
    lcp_array_.resize(size);
    vector<size_t> index_to_suffix(size);
    for (size_t i = 0; i < size; ++i) {
        index_to_suffix[suffix_array_[i]] = i;
    }
    size_t prev_lcp = 0;
    for (int current = 0; current < size; ++current) {
        if (index_to_suffix[current] == size - 1) {
            prev_lcp = 0;
        } else {
            if (prev_lcp > 0) {
                --prev_lcp;
            }
            size_t next = suffix_array_[index_to_suffix[current] + 1];
            while (current + prev_lcp < size
                   && next + prev_lcp < size
                   && string_[current + prev_lcp] == string_[next + prev_lcp]) {
                ++prev_lcp;
            }
        }
        lcp_array_[index_to_suffix[current]] = prev_lcp;
    }
}


class TwoSuffixTree {
public:
    TwoSuffixTree(const string &first, const string &second) {
        BuildFromTwoStrings(first, second);
    }

    void BuildFromTwoStrings(const string &first, const string &second);

    void PrintLexOrder(size_t node_id);

    size_t GetNodeSize() const {
        return trie_.size();
    }

private:
    struct Node {
        size_t parent;
        bool type;
        size_t left;
        size_t right;

        Node(size_t a, bool b, size_t c, size_t d) : parent(a), type(b), left(c), right(d) {}
    };

    void PrintLexOrderRec(size_t node_id);

    SuffixArrayWithLCP suf_;
    vector<Node> trie_;
    vector<vector<size_t>> children_;
    vector<size_t> lex_index_;
    size_t free_lex_ = 0;
};

void TwoSuffixTree::BuildFromTwoStrings(const string &first, const string &second) {
    string string = first + second;
    suf_ = SuffixArrayWithLCP(string);
    trie_.emplace_back(0, 0, 0, 0);
    size_t depth = 0;
    size_t current_node = 0;
    for (size_t i = 0; i < string.size(); ++i) {
        size_t lcp = i > 0 ? suf_.GetLCP(i - 1) : 0;
        size_t last_node;
        while (depth > lcp) {
            depth -= trie_[current_node].right - trie_[current_node].left;
            last_node = current_node;
            current_node = trie_[current_node].parent;
        }
        if (depth != lcp) {
            bool type = trie_[last_node].type;
            size_t left = trie_[last_node].left;
            size_t right = trie_[last_node].right;
            size_t middle = left - depth + lcp;

            trie_.emplace_back(current_node, type, left, middle);
            current_node = trie_.size() - 1;
            depth += trie_[current_node].right - trie_[current_node].left;

            trie_[last_node] = Node(trie_.size() - 1, type, middle, right);
        }
        size_t val = suf_[i] + lcp;
        bool type = suf_[i] >= first.size();
        size_t left = val - first.size() * type;
        size_t right = type ? second.size() : first.size();

        trie_.emplace_back(current_node, type, left, right);
        current_node = trie_.size() - 1;
        depth += trie_[current_node].right - trie_[current_node].left;
    }
    children_.resize(trie_.size(), vector<size_t>());
    for (size_t i = 1; i < trie_.size(); ++i) {
        children_[trie_[i].parent].push_back(i);
    }
}

void TwoSuffixTree::PrintLexOrderRec(size_t node_id) {
    if (node_id != 0 && lex_index_[node_id] == 0) {
        lex_index_[node_id] = free_lex_++;
    }
    if (node_id != 0) {
        Node node = trie_[node_id];
        cout << lex_index_[node.parent] << ' ' << node.type << ' ' << node.left << ' ' << node.right << '\n';
    }
    for (auto child_id : children_[node_id]) {
        PrintLexOrderRec(child_id);
    }
}

void TwoSuffixTree::PrintLexOrder(size_t node_id) {
    lex_index_.resize(trie_.size(), 0);
    free_lex_ = 1;
    PrintLexOrderRec(node_id);
}


class Solver {
public:
    void operator()() {
        string first, second;
        cin >> first >> second;
        TwoSuffixTree trie(first, second);
        cout << trie.GetNodeSize() << '\n';
        trie.PrintLexOrder(0);
    }
};

int main() {
    Solver solver;
    solver();

    return 0;
}