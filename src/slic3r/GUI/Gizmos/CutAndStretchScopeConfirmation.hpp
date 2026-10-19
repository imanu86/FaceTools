// Copyright (C) 2026 Emmanuele Semprini
// SPDX-License-Identifier: AGPL-3.0-only
// Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

#ifndef slic3r_CutAndStretchScopeConfirmation_hpp_
#define slic3r_CutAndStretchScopeConfirmation_hpp_

#include "libslic3r/CutAndStretch.hpp"

#include <memory>

namespace Slic3r::GUI {

// Consent authorizes one immutable prepared region, not an arbitrary later cut.
// A weak identity cannot keep old geometry alive or authorize an allocation that
// happens to reuse an old address. Source/pose/guide changes invalidate the plan
// in the real gizmo, whose completion path also rechecks every captured value.
class CutAndStretchScopeConfirmation {
public:
    void clear() noexcept { m_plan.reset(); }

    bool confirm(const std::shared_ptr<const CutAndStretchPlan> &plan,
                 CutAndStretchScope expected_scope) noexcept
    {
        clear();
        if (!plan || !*plan || plan->scope != expected_scope)
            return false;
        m_plan = plan;
        return true;
    }

    bool matches(const std::shared_ptr<const CutAndStretchPlan> &plan,
                 CutAndStretchScope expected_scope) const noexcept
    {
        return plan && bool(*plan) && plan->scope == expected_scope &&
               m_plan.lock() == plan;
    }

private:
    std::weak_ptr<const CutAndStretchPlan> m_plan;
};

} // namespace Slic3r::GUI

#endif // slic3r_CutAndStretchScopeConfirmation_hpp_
