#include <iostream>
#include <vector>
#include <string>

using std::string;
using std::vector;
using std::max;
using std::cout;
using std::cin;


class SuffixArrayWithLCP {
public:
    SuffixArrayWithLCP() = default;

    explicit SuffixArrayWithLCP(string line);

    size_t GetSize() const;

    size_t GetSuffix(size_t pos) const;

    size_t GetLCP(size_t pos) const;

private:
    void BuildSuffix();

    void BuildLCP();

    static const size_t alphabet_size = 256;
    string string_;
    vector<size_t> suffix_array_;
    vector<size_t> lcp_array_;
};

const size_t SuffixArrayWithLCP::alphabet_size;

SuffixArrayWithLCP::SuffixArrayWithLCP(string line) : string_(std::move(line)) {
    BuildSuffix();
    BuildLCP();
}

size_t SuffixArrayWithLCP::GetSize() const {
    return suffix_array_.size();
}

inline size_t SuffixArrayWithLCP::GetSuffix(size_t pos) const {
    return suffix_array_[pos];
}

inline size_t SuffixArrayWithLCP::GetLCP(size_t pos) const {
    return lcp_array_[pos];
}

void SuffixArrayWithLCP::BuildSuffix() {
    size_t size = string_.length();
    suffix_array_.assign(size, 0);

    vector<size_t> perm(size);
    vector<size_t> sorting_array(alphabet_size);
    vector<size_t> equivalence_classes(size);

    for (size_t i = 0; i < size; ++i) {
        sorting_array[string_[i]]++;
    }

    for (size_t i = 1; i < alphabet_size; ++i) {
        sorting_array[i] += sorting_array[i - 1];
    }

    for (size_t i = 0; i < size; ++i) {
        perm[--sorting_array[string_[i]]] = i;
    }

    equivalence_classes[perm[0]] = 0;
    size_t classes = 1;

    for (size_t i = 1; i < size; ++i) {
        if (string_[perm[i]] != string_[perm[i - 1]]) {
            ++classes;
        }
        equivalence_classes[perm[i]] = classes - 1;
    }

    vector<int> perm_tmp(size);
    vector<size_t> equivalence_classes_tmp(size);

    for (size_t i = 0; (1ULL << i) < size; ++i) {
        for (size_t j = 0; j < size; ++j) {
            perm_tmp[j] = perm[j] - (1ULL << i);
            if (perm_tmp[j] < 0) {
                perm_tmp[j] += size;
            }
        }

        sorting_array.assign(classes, 0);
        for (size_t j = 0; j < size; ++j) {
            sorting_array[equivalence_classes[perm_tmp[j]]]++;
        }

        for (size_t j = 1; j < classes; ++j) {
            sorting_array[j] += sorting_array[j - 1];
        }

        for (int j = size - 1; j >= 0; --j) {
            perm[--sorting_array[equivalence_classes[perm_tmp[j]]]] = perm_tmp[j];
        }

        equivalence_classes_tmp[perm[0]] = 0;
        classes = 1;

        for (size_t j = 1; j < size; ++j) {
            size_t mid1 = (perm[j] + (1ULL << i)) % size;
            size_t mid2 = (perm[j - 1] + (1ULL << i)) % size;
            if (equivalence_classes[perm[j]] != equivalence_classes[perm[j - 1]]
                || equivalence_classes[mid1] != equivalence_classes[mid2]) {
                classes++;
            }

            equivalence_classes_tmp[perm[j]] = classes - 1;
        }

        copy(equivalence_classes_tmp.begin(), equivalence_classes_tmp.begin() + size, equivalence_classes.begin());
    }

    copy(perm.begin(), perm.begin() + size, suffix_array_.begin());
}

void SuffixArrayWithLCP::BuildLCP() {
    size_t size = string_.length();
    lcp_array_.assign(size, 0);
    vector<size_t> pos(size);

    for (size_t i = 0; i < size; ++i) {
        pos[suffix_array_[i]] = i;
    }

    int k = 0;
    for (int i = 0; i < size; ++i) {
        if (k > 0) {
            k--;
        }

        if (pos[i] == size - 1) {
            lcp_array_[size - 1] = -1;
            k = 0;
        } else {
            int j = suffix_array_[pos[i] + 1];
            while (max(i + k, j + k) < size && string_[i + k] == string_[j + k]) {
                k++;
            }

            lcp_array_[pos[i]] = k;
        }
    }
}

size_t CountDifferentSubstrings(const string &line) {
    SuffixArrayWithLCP suffix_array_with_lcp(line);
    size_t ans = 0;
    size_t size = suffix_array_with_lcp.GetSize() - 1;

    for (size_t i = 0; i < size; ++i) {
        ans += size - suffix_array_with_lcp.GetSuffix(i) - suffix_array_with_lcp.GetLCP(i);
    }

    ans += size - suffix_array_with_lcp.GetSuffix(size);

    return ans;
}


int main() {
    string line;
    cin >> line;
    line += '$';

    cout << CountDifferentSubstrings(line);

    return 0;
}