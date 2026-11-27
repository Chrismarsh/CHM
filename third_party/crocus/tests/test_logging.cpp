#include "crocus_log.hpp"

#include <future>
#include <string>
#include <vector>

extern "C" void native_warning();

namespace
{
struct Capture
{
    std::vector<std::string>& warnings;
};

bool capture(crocus::LogLevel level, std::string_view text, const void* context)
{
    if (level == crocus::LogLevel::warning)
        static_cast<const Capture*>(context)->warnings.emplace_back(text);
    return true;
}
} // namespace

int main()
{
    std::vector<std::string> main_messages, worker_messages, nested_messages;
    const Capture main_log{main_messages}, worker_log{worker_messages}, nested_log{nested_messages};
    const crocus::LogScope main_scope(capture, &main_log);
    std::promise<void> ready, release;
    auto worker = std::async(std::launch::async,
                             [&]
                             {
                                 const crocus::LogScope worker_scope(capture, &worker_log);
                                 ready.set_value();
                                 release.get_future().wait();
                                 native_warning();
                             });
    ready.get_future().wait();
    native_warning(); // The worker has an active sink on its own thread.
    {
        const crocus::LogScope nested(capture, &nested_log);
        native_warning();
    }
    native_warning();
    release.set_value();
    worker.get();
    const std::string expected = "test native warning";
    return main_messages == std::vector<std::string>{expected, expected} &&
                   worker_messages == std::vector<std::string>{expected} && nested_messages == worker_messages
               ? 0
               : 1;
}
