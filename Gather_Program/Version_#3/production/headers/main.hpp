bool loop_simple_menu(int option_sk);
std::string get_value(const std::string msg);
void process_args(int argc, char* argv[]);
void inserts_ordered(int value, std::list<int> &lst_actns);
void exeute_command(std::list<int> &lst_actns, bool exe, std::string &vf, std::string &ti);
