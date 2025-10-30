# Security Policy

## Reporting a Vulnerability

Please do not open public issues for security vulnerabilities.

- Preferred: open a private advisory via GitHub Security Advisories ("Report a vulnerability" on the repo Security tab), or
- Email: <a.alzweidi@student.reading.ac.uk>

We aim to acknowledge reports within 72 hours and provide an initial assessment within 7 days. If the issue affects downstream users materially, we will coordinate a responsible disclosure and publish a patched release along with guidance.

## Scope

Vulnerabilities affecting confidentiality, integrity, or availability of the application or its users, including but not limited to:

- XSS/CSRF and injection attacks
- Supply-chain concerns (malicious dependencies)
- Leaks of secrets/tokens
- Logic bugs leading to data corruption or privilege misuse

## Out of scope (examples)

- DoS through unrealistic inputs that do not reflect supported environments
- Issues requiring privileged/local access without a clear privilege escalation
- Non-security bugs (please file a regular issue)

## Disclosure Process

1. Receive and confirm the report (triage).
2. Develop and validate a fix; prepare regression tests where feasible.
3. Coordinate a release and public disclosure notes.
4. Credit reporters who wish to be acknowledged.
