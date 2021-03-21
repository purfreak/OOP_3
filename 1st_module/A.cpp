#include <iostream>
#include <vector>

using std::string;
using std::vector;
using std::cout;
using std::cin;
using std::endl;


class PatternEntriesFinder {
public:
    explicit PatternEntriesFinder(string &line);

    void GetPatternEntries();

private:
    void CalculatePrefixFunction(const string &line);

    vector<size_t> prefix_function_array_;
    string pattern_;
};

PatternEntriesFinder::PatternEntriesFinder(string &line) : pattern_(line) {}

void PatternEntriesFinder::CalculatePrefixFunction(const string &line) {
    prefix_function_array_.resize(line.length(), 0);

    for (size_t i = 1; i < line.length(); ++i) {
        size_t last_value = prefix_function_array_[i - 1];

        while (last_value > 0 && line[i] != line[last_value]) {
            last_value = prefix_function_array_[last_value - 1];
        }

        if (line[i] == line[last_value]) {
            prefix_function_array_[i] = last_value + 1;
        }
    }
}

void PatternEntriesFinder::GetPatternEntries() {
    vector<size_t> pattern_entries;

    pattern_ += '#';
    CalculatePrefixFunction(pattern_);

    size_t last_value = 0;
    int i = 0;
    char next_symbol;

    while (cin >> next_symbol) {
        size_t curr_value = last_value;

        while (curr_value > 0 && next_symbol != pattern_[curr_value]) {
            curr_value = prefix_function_array_[curr_value - 1];
        }

        if (next_symbol == pattern_[curr_value]) {
            last_value = curr_value + 1;
        } else {
            last_value = 0;
        }

        if (last_value == pattern_.length() - 1) {
            pattern_entries.push_back(i - pattern_.length() + 2);
        }

        ++i;
    }

    for (auto entry : pattern_entries) {
        cout << entry << " ";
    }

    cout << endl;

}


int main() {
    string pattern;
    cin >> pattern;

    PatternEntriesFinder pattern_entries_finder(pattern);
    pattern_entries_finder.GetPatternEntries();

    return 0;
}
