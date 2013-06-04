---
title: CMS Offline Software
layout: default
related:
 - { name: Home, link: index.html }
 - { name: Project, link: https://github.com/cms-sw/cmssw }
 - { name: Topic Collector, link: https://cern.ch/cmsgit/cmsgit}
 - { name: Feedback, link: https://github.com/cms-sw/cmssw/issues/new }
---

## Sparse checkout does not work.

Apparently some university deployed a non working `git` 1.7.4 client. This
results in sparse checkout misbehavior. Using 1.7.4.1 or later seems to fix the
issue.
