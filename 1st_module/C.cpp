#include <memory>
#include <vector>
#include <iostream>

using std::cin;
using std::cout;
using std::string;
using std::shared_ptr;
using std::weak_ptr;
using std::vector;
using std::pair;
using std::make_shared;

const size_t alphabet = 26;

struct Node {
public:
    Node(const shared_ptr<Node> &parent, char symbol);

    shared_ptr<Node> GetParent() const;

    weak_ptr<Node> parent;
    shared_ptr<Node> suff_link;
    shared_ptr<Node> up_link;
    vector<shared_ptr<Node>> vertices;
    vector<size_t> patterns;
    bool terminal;
    char edge_symbol;
};

class AhoCorasic {
public:
    explicit AhoCorasic(string &mask);

    vector<size_t> GetEntries(const string &text);

private:
    void AddPattern(size_t index);

    shared_ptr<Node> GetSonBySymbol(const shared_ptr<Node> &node, char symbol);

    shared_ptr<Node> GetLink(const shared_ptr<Node> &node, char symbol);

    shared_ptr<Node> GetSuffLink(shared_ptr<Node> node);

    shared_ptr<Node> GetUpLink(const shared_ptr<Node> &node);

    shared_ptr<Node> root_;
    vector<pair<size_t, size_t>> patterns_;
    string mask_;
};

Node::Node(const shared_ptr<Node> &parent, char symbol)
        : parent(parent), edge_symbol(symbol), vertices(alphabet),
          terminal(false) {}

shared_ptr<Node> Node::GetParent() const {
    return parent.lock();
}

AhoCorasic::AhoCorasic(string &mask) : root_(make_shared<Node>(nullptr, EOF)),
                                       mask_(mask) {
    int pattern_begin = 0;
    int pattern_end = 0;

    while (pattern_end < mask.size()) {
        pattern_end = mask.find('?', pattern_begin);

        if (pattern_end == -1) {
            pattern_end = mask.size();
        }

        if (pattern_end - pattern_begin > 0) {
            patterns_.emplace_back(pattern_begin, pattern_end);
        }

        pattern_begin = pattern_end + 1;
    }

    for (int i = 0; i < patterns_.size(); ++i) {
        AddPattern(i);
    }
}

vector<size_t> AhoCorasic::GetEntries(const string &text) {
    auto current_node = root_;
    auto entries = vector<int>(text.size());

    for (int i = 0; i < text.size(); ++i) {
        char current_symbol = text[i];
        current_node = GetLink(current_node, current_symbol);

        auto up_link = current_node;

        while (up_link != root_) {
            if (up_link->terminal) {
                for (auto &pattern : up_link->patterns) {
                    if (i + 1 >= patterns_[pattern].second) {
                        ++entries[i + 1 - patterns_[pattern].second];
                    }
                }
            }
            up_link = GetUpLink(up_link);
        }
    }

    vector<size_t> answer;

    for (int i = 0; i < entries.size(); ++i) {
        if (entries[i] == patterns_.size()) {
            if (i + mask_.size() <= text.size()) {
                answer.push_back(i);
            }
        }
    }

    return answer;
}

void AhoCorasic::AddPattern(size_t index) {
    auto pattern_position = patterns_[index];
    auto current_node = root_;

    for (int i = pattern_position.first; i < pattern_position.second; ++i) {
        char symbol = mask_[i];
        auto node = GetSonBySymbol(current_node, symbol);

        if (node == nullptr) {
            current_node->vertices[symbol - 'a'] =
                    make_shared<Node>(current_node, symbol);
        }

        current_node = current_node->vertices[symbol - 'a'];
    }

    current_node->terminal = true;
    current_node->patterns.push_back(index);
}

shared_ptr<Node> AhoCorasic::GetSonBySymbol(const shared_ptr<Node> &node,
                                            char symbol) {
    if (symbol == EOF && node == root_) {
        return root_;
    }

    return node->vertices[symbol - 'a'];
}

shared_ptr<Node> AhoCorasic::GetLink(const shared_ptr<Node> &node,
                                     char symbol) {
    if (symbol == EOF && node == root_) {
        return root_;
    }

    if (GetSonBySymbol(node, symbol) != nullptr) {
        return GetSonBySymbol(node, symbol);
    }

    if (node == root_) {
        node->vertices[symbol - 'a'] = node;
    } else {
        node->vertices[symbol - 'a'] =
                GetLink(GetSuffLink(node), symbol);
    }

    return GetSonBySymbol(node, symbol);
}

shared_ptr<Node> AhoCorasic::GetSuffLink(shared_ptr<Node> node) {
    if (node->suff_link == nullptr) {
        auto parent = node->GetParent();

        if (parent == nullptr) {
            return node;
        }

        if (node == root_ || parent == root_) {
            node->suff_link = node == root_ ? node : parent;
        } else {
            node->suff_link =
                    GetLink(GetSuffLink(parent), node->edge_symbol);
        }
    }

    return node->suff_link;
}

shared_ptr<Node> AhoCorasic::GetUpLink(const shared_ptr<Node> &node) {
    if (node->up_link == nullptr) {
        if (GetSuffLink(node)->terminal) {
            node->up_link = GetSuffLink(node);
        } else if (GetSuffLink(node) == root_) {
            node->up_link = root_;
        } else {
            node->up_link = GetUpLink(GetSuffLink(node));
        }
    }

    return node->up_link;
}

class Solver {
public:
    void operator()() {
        string text;
        string pattern;

        cin >> pattern >> text;

        AhoCorasic aho_corasic(pattern);

        for (auto entry: aho_corasic.GetEntries(text)) {
            cout << entry << " ";
        }
    }
};

int main() {
    Solver solver;
    solver();
}
