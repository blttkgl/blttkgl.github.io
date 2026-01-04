---
title: 'Automation beats perfection'
date: 2026-01-06
permalink: /posts/2026/01/automation_beats_perfection/
tags:
  - OpenFOAM
  - CFD
  - Automation
  - Simulation Workflow
  - Reproducibility
---

One of the most counterintuitive lessons I learned over time is that **a slightly imperfect CFD setup that runs reliably** will outperform a "perfect" setup that needs constant babysitting.

Early in my career I spent an enormous amount of time to make individual cases beautiful: carefully thought out solver settings, well designed perfect meshes, carefully adjusted time step controls. The result was often a single impressible run, and weeks of pain when anything changed in the setup.

More often than not, **CFD is not an one-off exercise.**

Instead, you rerun cases because:

- Geometry changes
- Boundary conditions change
- Models change
- Meshes change
- Someone asks "what if?"

If any of these changes require manual intervention, your workflow becomes fragile very quickly. 

Automation helps you **contain the complexity.**

Many clever "CFD tricks" utilize:

- Hand picked time stepping strategy
- Case specific tweaks
- One-off scripts that does the trick but very hard to maintain and expand

They work because everything else stays exactly the same. As soon as you:

- Refine the mesh
- Change input parameters (e.g., valve or injection timings)
- Switch models (e.g., turbulence, injection)

those tricks stop working. 

A robust, automated workflow **survives change.**

## Reproducibility is an engineering requirement, not a luxury

In CFD setups it is very easy to lose track of **what actually produced the result**.

Without automation:

- Dictionaries drift
- Scripts get edited manually
- "final_final_final_v4" cases appear

Automation forces discipline:

- Fully connected case setups where changes in user input is reflected everywhere
- Consistent execution
- Trecable output

If you cannot rerun your case six months later and get the same result, you don't have a workflow, you have an one-shot case.

## "Good enough, repeatable" beats "perfect, once"

In industrial CFD cases especially, the value comes from trends, sensitivities, and comparisons.

 An acceptable case setup that runs every time, responds consistently to changes, and finished within a predictable time is often more valuable than a fragile high-fidelity run that fails half the time, or runs for much longer than you can afford.

## What automation looks like in practice

For me, this eventually meant:

- Scripted meshing and case setup (which I will talk about in a later post)
- Standardized directory structures
- Parametrized dictionaries (again, another post 😊)
- Effective and well thought out run-time processing that generates all the useful data you need
- Minimal manual editing

## If–else conditionals and macros are your friend

One of the most underused features in OpenFOAM case setups is the ability to use **macros and conditional logic** directly in dictionaries.

Early on, I treated dictionaries as static configuration files. Over time, I realized they can (and should) behave more like **parameterized inputs**.

Using:
- `#include`, `#includeIfPresent`
- `#calc`
- `#if`, `#else`, `#endif`, `#ifEq`
- macros

you can build case setups that **adapt automatically** instead of being manually duplicated and edited.

Typical examples where this helps enormously:
- switching models on and off without duplicating dictionaries
- changing time-step controls based on solver or physics
- reusing the same case structure for multiple geometries
- avoiding copy–paste errors across similar cases

Instead of maintaining:
> `case_LES/`, `case_RANS/`, `case_cold/`, `case_reactive/`

you maintain **one case**, driven by a small set of high-level parameters.

This does two important things:
- it reduces human error
- it makes your assumptions explicit and traceable

Macros and conditionals won’t make a bad setup good but they make a good setup **scalable, reproducible, and robust**.

If you find yourself copying dictionaries and editing the same lines over and over, that’s usually a sign that logic belongs in the case, and not in your memory.

I strongly recommend checking out the ["Basic input/output file format" section of the OpenFOAM User Guide](https://doc.cfd.direct/openfoam/user-guide-v13/basic-file-format).

The biggest mental shift is this:

> **Stop optimizing the case.**  
> **Start optimizing the process.**

Because in the long run, the workflow that survives change is the one that produces insight. 

Thanks for reading,  
— Bulut
