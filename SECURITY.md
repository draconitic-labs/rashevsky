# Security

Report vulnerabilities privately through [GitHub private vulnerability reporting](https://github.com/draconitic-labs/rashevsky/security/advisories/new), or at [support@draconitic.com](mailto:support@draconitic.com).

Use the latest release. Fixes go there, and there are no maintained older branches.

The desktop application downloads published models only when asked, over TLS, from the addresses built into it, and refuses a file whose SHA-256 digest differs from the one built into it. SBML files you open yourself are parsed in-process, so open models only from sources you trust.
