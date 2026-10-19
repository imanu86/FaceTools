// Copyright (C) 2026 Emmanuele Semprini
// SPDX-License-Identifier: AGPL-3.0-only
// Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

#include <catch2/catch.hpp>
#include "slic3r/GUI/Jobs/Job.hpp"
#include <wx/app.h>
#include <wx/init.h>
#include <wx/thread.h>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {
// Dedicated test process only. This is the native wx console runtime/event
// queue, NOT GUI_App: no windows, profiles, printers, network or plugins.
// Do not call Job::start() or emit terminal progress (busy-cursor GUI paths).
class WxConsoleRuntime {
    int argc = 1;
    char name[32] = "job-ui-dispatch-tests";
    char *argv[2] = {name, nullptr};
public:
    const std::thread::id main_thread = std::this_thread::get_id();
    WxConsoleRuntime()
    {
        if (wxAppConsole::GetInstance())
            throw std::runtime_error("Job dispatch tests require a dedicated process without an existing wx app");
        wxAppConsole::SetInstance(new wxAppConsole);
        if (!wxEntryStart(argc, argv))
            throw std::runtime_error("Cannot initialize native wx console runtime");
        if (!wxAppConsole::GetInstance()->CallOnInit()) {
            wxEntryCleanup();
            throw std::runtime_error("Cannot initialize wx console app");
        }
    }
    ~WxConsoleRuntime()
    {
        wxAppConsole::GetInstance()->ProcessPendingEvents();
        wxAppConsole::GetInstance()->OnExit();
        wxEntryCleanup();
    }
    void drain()
    {
        REQUIRE(std::this_thread::get_id() == main_thread);
        REQUIRE(wxIsMainThread());
        wxAppConsole::GetInstance()->ProcessPendingEvents();
    }
};
WxConsoleRuntime &runtime()
{
    static WxConsoleRuntime instance;
    return instance;
}

struct Call {
    std::string kind, text, description, extra;
    int value = 0;
    std::thread::id thread;
    bool wx_main = false;
};
class CallLog {
    mutable std::mutex mutex;
    std::vector<Call> calls;
public:
    void add(std::string kind, int value = 0, std::string text = {}, std::string description = {}, std::string extra = {})
    {
        std::lock_guard<std::mutex> lock(mutex);
        calls.push_back({std::move(kind), std::move(text), std::move(description), std::move(extra), value,
                         std::this_thread::get_id(), wxIsMainThread()});
    }
    std::vector<Call> snapshot() const
    {
        std::lock_guard<std::mutex> lock(mutex);
        return calls;
    }
};
class FakeProgress : public Slic3r::ProgressIndicator {
    std::shared_ptr<CallLog> log;
public:
    explicit FakeProgress(std::shared_ptr<CallLog> target) : log(std::move(target)) {}
    void clear_percent() override { log->add("clear"); }
    void show_error_info(wxString msg, int code, wxString description, wxString extra) override
    { log->add("error", code, msg.utf8_string(), description.utf8_string(), extra.utf8_string()); }
    void set_range(int value) override { log->add("range", value); }
    void set_cancel_callback(CancelFn = CancelFn()) override { log->add("cancel_callback"); }
    void set_progress(int value) override { log->add("progress", value); }
    void set_status_text(const char *text) override { log->add("status", 0, text ? text : ""); }
    int get_range() const override { log->add("get_range"); return 100; }
};
class ProbeJob : public Slic3r::GUI::Job {
public:
    using Job::Job;
    void emit_error(wxString message, int code, wxString description = {}, wxString extra = {})
    { show_error_info(message, code, description, extra); }
    void emit_status(int value, const wxString &message) { update_status(value, message); }
    void mark_finalized() { Job::finalize(); }
    bool is_print_job() const override { return true; }
protected:
    void process() override {} // never launched; the test owns/joins its std::thread
};
struct Fixture {
    WxConsoleRuntime &wx = runtime();
    std::shared_ptr<CallLog> log = std::make_shared<CallLog>();
    std::shared_ptr<FakeProgress> progress = std::make_shared<FakeProgress>(log);
    std::unique_ptr<ProbeJob> job = std::make_unique<ProbeJob>(progress);
};
void all_on_main(const std::vector<Call> &calls)
{
    for (const auto &call : calls) {
        CHECK(call.thread == runtime().main_thread);
        CHECK(call.wx_main);
    }
}
}

TEST_CASE("Job errors are deferred to the native wx main thread with owned complete payloads", "[job_ui_dispatch]")
{
    Fixture f;
    const std::string message = "Synthetic upload error: caffè 漢字";
    const std::string description = "Synthetic description: λ";
    const std::string extra = "Synthetic extra: details only";
    std::thread worker([&] {
        wxString msg = wxString::FromUTF8(message), desc = wxString::FromUTF8(description), info = wxString::FromUTF8(extra);
        f.job->emit_error(msg, -7301, desc, info);
        msg.assign("changed message"); desc.clear(); info.assign("changed extra");
    });
    worker.join();
    CHECK(f.log->snapshot().empty()); // nonfatal: old Job exposes its worker-thread call
    f.wx.drain();
    const auto calls = f.log->snapshot();
    REQUIRE(calls.size() == 1);
    CHECK(calls[0].kind == "error");
    CHECK(calls[0].text == message);
    CHECK(calls[0].value == -7301);
    CHECK(calls[0].description == description);
    CHECK(calls[0].extra == extra);
    all_on_main(calls);
}

TEST_CASE("Job status error and status retain FIFO ordering in the same event target", "[job_ui_dispatch]")
{
    Fixture f;
    std::thread worker([&] {
        f.job->emit_status(10, "before");
        f.job->emit_error("failure", 42, "description", "extra");
        f.job->emit_status(20, "after");
    });
    worker.join();
    CHECK(f.log->snapshot().empty());
    f.wx.drain();
    const auto calls = f.log->snapshot();
    REQUIRE(calls.size() == 5);
    CHECK(calls[0].kind == "status"); CHECK(calls[0].text == "before");
    CHECK(calls[1].kind == "progress"); CHECK(calls[1].value == 10);
    CHECK(calls[2].kind == "error"); CHECK(calls[2].value == 42);
    CHECK(calls[3].kind == "status"); CHECK(calls[3].text == "after");
    CHECK(calls[4].kind == "progress"); CHECK(calls[4].value == 20);
    all_on_main(calls);
}

TEST_CASE("Job cancellation suppresses both queued errors and errors emitted after cancellation", "[job_ui_dispatch]")
{
    Fixture f;
    std::thread before_cancel([&] { f.job->emit_error("queued before cancel", 1); });
    before_cancel.join();
    f.job->cancel();
    CHECK(f.log->snapshot().empty());
    f.wx.drain();
    CHECK(f.log->snapshot().empty());
    std::thread after_cancel([&] { f.job->emit_error("late canceled error", 2); });
    after_cancel.join();
    f.wx.drain();
    CHECK(f.log->snapshot().empty());
}

TEST_CASE("Job finalization alone does not hide queued or later print upload errors", "[job_ui_dispatch]")
{
    Fixture f;
    std::thread first([&] { f.job->emit_error("error queued before finalize", 31); });
    first.join();
    f.job->mark_finalized(); // no status=100/busy cursor: test the real finalized flag directly
    REQUIRE(f.job->is_finalized());
    CHECK(f.log->snapshot().empty());
    f.wx.drain();
    auto calls = f.log->snapshot();
    REQUIRE(calls.size() == 1);
    CHECK(calls[0].text == "error queued before finalize");
    std::thread second([&] { f.job->emit_error("failure after waiting at 100 percent", 32); });
    second.join();
    CHECK(f.log->snapshot().size() == 1);
    f.wx.drain();
    calls = f.log->snapshot();
    REQUIRE(calls.size() == 2);
    CHECK(calls[1].value == 32);
    CHECK(calls[1].text == "failure after waiting at 100 percent");
    all_on_main(calls);
}

TEST_CASE("Job destruction discards pending native events before releasing its progress owner", "[job_ui_dispatch]")
{
    Fixture f;
    std::weak_ptr<FakeProgress> weak = f.progress;
    std::thread worker([&] {
        f.job->emit_error("must never reach destroyed progress", 71);
        f.job->emit_status(7, "also discarded");
    });
    worker.join(); // required owner contract: NEVER destroy while producers still run
    CHECK(f.log->snapshot().empty());
    f.progress.reset();
    CHECK_FALSE(weak.expired()); // Job owns the remaining progress reference
    f.job.reset();
    CHECK(weak.expired());
    f.wx.drain();
    CHECK(f.log->snapshot().empty());
}

TEST_CASE("Job retry uses a fresh target and cannot deliver the retired job error to shared progress", "[job_ui_dispatch]")
{
    Fixture f;
    std::thread old_worker([&] { f.job->emit_error("old retired job", 81); });
    old_worker.join();
    f.job->cancel();
    f.job.reset();
    f.job = std::make_unique<ProbeJob>(f.progress); // same progress, new Job as SelectMachine retry does
    std::thread new_worker([&] { f.job->emit_error("new active job", 82, "new description", "new extra"); });
    new_worker.join();
    CHECK(f.log->snapshot().empty());
    f.wx.drain();
    const auto calls = f.log->snapshot();
    REQUIRE(calls.size() == 1);
    CHECK(calls[0].text == "new active job");
    CHECK(calls[0].value == 82);
    CHECK(calls[0].description == "new description");
    CHECK(calls[0].extra == "new extra");
    all_on_main(calls);
}

TEST_CASE("Job error events remain target specific when another live job is canceled", "[job_ui_dispatch]")
{
    Fixture first, second;
    std::thread worker([&] {
        first.job->emit_error("canceled target", 91);
        for (int i = 0; i < 8; ++i)
            second.job->emit_error(wxString::Format("active target %d", i), 1000 + i);
    });
    worker.join();
    first.job->cancel();
    CHECK(first.log->snapshot().empty());
    CHECK(second.log->snapshot().empty());
    first.wx.drain();
    CHECK(first.log->snapshot().empty());
    const auto calls = second.log->snapshot();
    REQUIRE(calls.size() == 8);
    for (int i = 0; i < 8; ++i) {
        CHECK(calls[i].value == 1000 + i);
        CHECK(calls[i].text == "active target " + std::to_string(i));
    }
    all_on_main(calls);
}
