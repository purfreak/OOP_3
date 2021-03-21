#include <algorithm>
#include <cmath>
#include <iostream>
#include <queue>
#include <vector>

using std::cin;
using std::cout;
using std::istream;
using std::max;
using std::queue;
using std::string;
using std::vector;

constexpr int ascii_shift_for_letters = 96;
constexpr int ascii_shift_for_digits = 48;
constexpr int max_possible_positions = 8192;

constexpr int left_border = 1;
constexpr int right_border = 8;

constexpr int KW1 = 3;
constexpr int KW2 = 3;

class Endspiel {
public:
    struct State {
        State();

        int KB_1;
        int KB_2;
        int QW_1;
        int QW_2;
        bool white_turn_to_play;

        bool IsCorrect();

        bool IsCheck();

        bool IsCheckmate();

        int CodeFromPosition();

        friend istream &operator>>(istream &is, State &state);
    };

    explicit Endspiel(State start_);

    static State PositionFromCode(int code);

    int CalculateMovesToCheckmate();

private:
    int start_index_;
    vector<bool> visited_ = vector<bool>(max_possible_positions, false);
    vector<int> moves_to_checkmate_ = vector<int>(max_possible_positions, 0);

    static void GetChildrenOrAncestors(int position_code, vector<int> &children,
                                       bool is_search_for_children);
};

Endspiel::Endspiel(State start_) : start_index_(start_.CodeFromPosition()) {}

Endspiel::State::State()
        : KB_1(0), KB_2(0), QW_1(0), QW_2(0), white_turn_to_play(true) {}

istream &operator>>(istream &is, Endspiel::State &state) {
    string QW, KB;
    is >> QW >> KB;
    state.QW_1 = QW[0] - ascii_shift_for_letters;
    state.QW_2 = QW[1] - ascii_shift_for_digits;
    state.KB_1 = KB[0] - ascii_shift_for_letters;
    state.KB_2 = KB[1] - ascii_shift_for_digits;
    return is;
}

bool Endspiel::State::IsCorrect() {
    if (QW_1 < left_border || QW_2 < left_border || KB_1 < left_border || KB_2 < left_border) {
        return false;
    }

    if (QW_1 > right_border || QW_2 > right_border || KB_1 > right_border || KB_2 > right_border) {
        return false;
    }

    // when KB is next to KW
    if (KB_2 <= (KW2 + 1) && KB_2 >= (KW2 - 1) && KB_1 <= (KW1 + 1) && KB_1 >= (KW1 - 1)) {
        return false;
    }

    if (QW_2 == KB_2 && QW_1 == KB_1) {
        return false;
    }

    if (QW_2 == KW2 && QW_1 == KW1) {
        return false;
    }

    if (white_turn_to_play) {
        if (IsCheck()) {
            return false;
        }
    }

    // when QW and KB are on the cells-neighbours
    if (max(abs(QW_2 - KB_2), abs(QW_1 - KB_1)) == 1) {
        // when KB made it a check
        if (white_turn_to_play) {
            return false;
        }
        // we don't want it to become a draw
        if (max(abs(QW_2 - KW2), abs(QW_1 - KW1)) != 1) {
            return false;
        }
    }

    return true;
}

bool Endspiel::State::IsCheck() {
    if (QW_1 == KB_1) {
        if (QW_1 != KW1) {
            return true;
        }
        return (KB_2 < KW2 && QW_2 < KW2) || (KB_2 > KW2 && QW_2 > KW2);
    }

    if (QW_2 == KB_2) {
        if (QW_2 != KW2) {
            return true;
        }
        return (KB_1 < KW1 && QW_1 < KW1) || (KB_1 > KW1 && QW_1 > KW1);
    }

    if (QW_1 - KB_1 == QW_2 - KB_2) {
        if (QW_1 - KW1 != QW_2 - KW2) {
            return true;
        }
        return ((QW_2 > KW2 && KB_2 > KW2) || (QW_2 < KW2 && KB_2 < KW2));
    }

    if (QW_1 - KB_1 == KB_2 - QW_2) {
        if (QW_1 - KW1 != KW2 - QW_2) {
            return true;
        }
        return ((QW_1 > KW1 && KB_1 > KW1) || (QW_1 < KW1 && KB_1 < KW1));
    }

    return false;
}

bool Endspiel::State::IsCheckmate() {
    if (white_turn_to_play || !IsCheck()) {
        return false;
    }
    vector<int> positions;
    GetChildrenOrAncestors(CodeFromPosition(), positions, true);
    return positions.empty();
}

void Endspiel::GetChildrenOrAncestors(int position_code, vector<int> &children,
                                      bool is_search_for_children) {
    State position = PositionFromCode(position_code);
    if ((position.white_turn_to_play && is_search_for_children) ||
        (!position.white_turn_to_play && !is_search_for_children)) {
        position.white_turn_to_play = !position.white_turn_to_play;
        for (int i = -7; i < 8; ++i) {
            if (i != 0) {
                vector<State> possible_positions = vector<State>(4, position);

                possible_positions[0].QW_1 += i;
                possible_positions[1].QW_2 += i;
                possible_positions[2].QW_1 += i;
                possible_positions[2].QW_2 += i;
                possible_positions[3].QW_1 += i;
                possible_positions[3].QW_2 -= i;

                for (auto pos : possible_positions) {
                    if (pos.IsCorrect()) {
                        double first_distance = hypot(abs(position.QW_1 - pos.QW_1),
                                                      abs(position.QW_2 - pos.QW_2));
                        double second_distance =
                                hypot(abs(position.QW_1 - KW1), abs(position.QW_2 - KW2)) +
                                hypot(abs(pos.QW_1 - KW1), abs(pos.QW_2 - KW2));
                        //check if the king isn't between QW's new position and QW's old position
                        if (fabs(first_distance - second_distance) > 1e-4) {
                            children.push_back(pos.CodeFromPosition());
                        }
                    }
                }
            }
        }
    } else {
        position.white_turn_to_play = !position.white_turn_to_play;
        for (int i = -1; i < 2; ++i) {
            for (int j = -1; j < 2; ++j) {
                if (i != 0 || j != 0) {
                    State current_position = position;
                    current_position.KB_1 += i;
                    current_position.KB_2 += j;
                    if (current_position.IsCorrect()) {
                        children.push_back(current_position.CodeFromPosition());
                    }
                }
            }
        }
    }
}

int Endspiel::State::CodeFromPosition() {
    int position_number =
            (KB_1 - 1) + ((KB_2 - 1) << 3) + ((QW_1 - 1) << 6) + ((QW_2 - 1) << 9);
    return white_turn_to_play ? position_number : position_number + 4096;
}

Endspiel::State Endspiel::PositionFromCode(int code) {
    Endspiel::State current_position;
    if (code >= 4096) {
        code -= 4096;
        current_position.white_turn_to_play = false;
    } else {
        current_position.white_turn_to_play = true;
    }

    current_position.QW_2 = (code - code % (1 << 9)) / (1 << 9) + 1;
    code = code % (1 << 9);
    current_position.QW_1 = (code - code % (1 << 6)) / (1 << 6) + 1;
    code = code % (1 << 6);
    current_position.KB_2 = (code - code % (1 << 3)) / (1 << 3) + 1;
    code = code % (1 << 3);
    current_position.KB_1 = code + 1;
    return current_position;
}

int Endspiel::CalculateMovesToCheckmate() {
    queue<int> positions_to_scan_indexes;
    for (int i = 0; i < max_possible_positions; ++i) {
        State current_position = PositionFromCode(i);
        if (current_position.IsCorrect() && current_position.IsCheckmate()) {
            positions_to_scan_indexes.push(i);
            moves_to_checkmate_[i] = 0;
            visited_[i] = true;
        }
    }
    while (!positions_to_scan_indexes.empty()) {
        int cur_position_index = positions_to_scan_indexes.front();
        positions_to_scan_indexes.pop();
        State current_position = PositionFromCode(cur_position_index);

        vector<int> ancestors_indexes;
        GetChildrenOrAncestors(cur_position_index, ancestors_indexes, false);

        if (!current_position.white_turn_to_play) {
            for (auto ancestor_index : ancestors_indexes) {
                if (!visited_[ancestor_index]) {
                    moves_to_checkmate_[ancestor_index] = moves_to_checkmate_[cur_position_index] + 1;
                    positions_to_scan_indexes.push(ancestor_index);
                    visited_[ancestor_index] = true;
                    if (ancestor_index == start_index_) {
                        return moves_to_checkmate_[start_index_];
                    }
                }
            }
        } else {
            for (auto ancestor_index : ancestors_indexes) {
                if (!visited_[ancestor_index]) {
                    int moves_to_checkmate = -1;
                    bool ready_to_be_pushed = true;

                    vector<int> children_indexes;
                    GetChildrenOrAncestors(ancestor_index, children_indexes, true);

                    for (auto child_index : children_indexes) {
                        if (visited_[child_index]) {
                            moves_to_checkmate =
                                    max(moves_to_checkmate, moves_to_checkmate_[child_index] + 1);
                        } else {
                            ready_to_be_pushed = false;
                            break;
                        }
                    }
                    if (ready_to_be_pushed) {
                        moves_to_checkmate_[ancestor_index] = moves_to_checkmate;
                        positions_to_scan_indexes.push(ancestor_index);
                        visited_[ancestor_index] = true;
                    }
                }
            }
        }
    }

    return moves_to_checkmate_[start_index_];
}

class Solver {
public:
    void operator()() {
        Endspiel::State state;
        cin >> state;
        Endspiel endspiel(state);
        cout << endspiel.CalculateMovesToCheckmate();
    }
};

int main() {
    Solver solver;
    solver();

    return 0;
}