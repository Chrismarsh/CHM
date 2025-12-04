#include <concepts>

template<typename T, typename Data>
concept GuaranteeImplementExecute = requires(T t, Data d) {
    { t.execute_impl(d) } -> std::same_as<void>;
};

template<class Derived, class Data>
    requires GuaranteeImplementExecute<Derived, Data>
class base_step {
public:
    void execute(Data& d) {
        static_cast<Derived*>(this)->execute_impl(d);
    }
};
