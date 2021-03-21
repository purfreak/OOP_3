#include <iostream>
#include <vector>

using std::cin;
using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::min;
using std::max;


vector<size_t> ZFunctionByLine(const string &line) {
    vector<size_t> z_function_array(line.length(), 0);

    size_t left = 0;
    size_t right = 0;
    for (size_t i = 1; i < line.length(); ++i) {
        size_t temp = static_cast<size_t>(max(0, min(static_cast<int>(right - i),
                                                     static_cast<int>(z_function_array[i - left]))));

        while (i + temp < line.size() && line[temp] == line[i + temp]) {
            temp++;
        }

        if (i + temp > right) {
            left = i;
            right = i + temp;
        }

        z_function_array[i] = temp;
    }

    return z_function_array;
}

string LineByZFunction(const vector<size_t> &z_function_array) {
    if (z_function_array.empty()) {
        return "";
    }

    vector<vector<size_t>> vec(z_function_array.size());

    for (size_t i = 1; i < z_function_array.size(); ++i) {
        if (!z_function_array[i]) {
            vec[z_function_array[i] + i - 1].push_back(z_function_array[i]);
        }
    }

    string result_line = "a";

    for (size_t i = 1; i < z_function_array.size(); ++i) {
        static size_t prefix_size = 0;
        static size_t prefix_move = 0;

        if (z_function_array[i] == 0 && prefix_size == 0) {
            size_t vec_count = 0;
            vector<bool> used('z' - 'a', false);
            used[result_line[vec_count] - 'a'] = true;

            while (vec_count < vec[i - 1].size()) {
                used[result_line[vec[i - 1][vec_count]] - 'a'] = true;
                vec_count++;
            }

            char next_char = 'a';
            for (auto u : used) {
                if (!u) {
                    result_line += next_char;
                    break;
                }
                next_char++;
            }

        }

        if (z_function_array[i] > prefix_size) {
            prefix_size = z_function_array[i];
            prefix_move = 0;
        }

        if (prefix_size > 0) {
            result_line += result_line[prefix_move];
            prefix_size--;
            prefix_move++;
        }
    }

    return result_line;
}

vector<size_t> PrefixFunctionByLine(const string &line) {
    vector<size_t> prefix_function_array(line.length(), 0);

    for (size_t i = 1; i < line.size(); ++i) {
        size_t curr_value = prefix_function_array[i - 1];

        while (curr_value > 0 && line[i] != line[curr_value]) {
            curr_value = prefix_function_array[curr_value - 1];
        }

        if (line[i] == line[curr_value]) {
            prefix_function_array[i] = curr_value + 1;
        }
    }

    return prefix_function_array;
}

string LineByPrefixFunction(const vector<size_t> &prefix_function_array) {
    if (prefix_function_array.empty()) {
        return "";
    }

    string result_line = "a";
    for (size_t i = 1; i < prefix_function_array.size(); ++i) {
        if (prefix_function_array[i] != 0) {
            result_line += result_line[prefix_function_array[i] - 1];
        } else {
            vector<bool> used('z' - 'a', false);
            used[0] = true;
            size_t line_count = prefix_function_array[i - 1];
            used[result_line[line_count] - 'a'] = true;

            while (line_count != 0) {
                used[result_line[line_count] - 'a'] = true;
                line_count = prefix_function_array[line_count - 1];
            }

            char next_char = 'a';
            for (auto u : used) {
                if (!u) {
                    result_line += next_char;
                    break;
                }
                next_char++;
            }
        }
    }

    return result_line;
}

vector<size_t> PrefixFunctionByZFunction(const vector<size_t> &z_function_array) {
    vector<size_t> prefix_function_array(z_function_array.size(), 0);

    for (size_t i = 1; i < z_function_array.size(); ++i) {
        for (size_t j = z_function_array[i] - 1; j > 0; --j) {
            if (prefix_function_array[i + j] > 0) {
                break;
            } else {
                prefix_function_array[i + j] = j + 1;
            }
        }
    }

    return prefix_function_array;
}

vector<size_t> ZFunctionByPrefixFunction(const vector<size_t> &prefix_function_array) {
    vector<size_t> z_function_array(prefix_function_array.size());

    for (size_t i = 1; i < prefix_function_array.size(); ++i) {
        if (prefix_function_array[i] > 0) {
            z_function_array[i - prefix_function_array[i] + 1] = prefix_function_array[i];
        }
    }

    z_function_array[0] = prefix_function_array.size();
    for (size_t i = 1; i < prefix_function_array.size();) {
        size_t z_position = i;

        if (z_function_array[i] > 0) {
            for (size_t z_move = 1; z_move < z_function_array[i]; z_move++) {
                if (z_function_array[i + z_move] > z_function_array[z_move]) {
                    break;
                }

                z_function_array[i + z_move] = std::min(z_function_array[z_move], z_function_array[i] - z_move);
                z_position = i + z_move;
            }
        }
        i = z_position + 1;
    }

    z_function_array[0] = 0;

    return z_function_array;
}


int main() {
    size_t symbol;
    vector<size_t> prefix_function_array;

    while (cin >> symbol) {
        prefix_function_array.push_back(symbol);
    }

    cout << LineByPrefixFunction(prefix_function_array) << endl;
}
