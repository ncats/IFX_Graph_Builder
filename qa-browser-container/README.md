# QA Browser container

The image needs read-only access to the internal `IFX_Registry` repository while
it is built.

For a local build, the base Compose configuration forwards your SSH agent.
Load an identity with access to `ncats/IFX_Registry` if needed, then build:

```bash
ssh-add ~/.ssh/id_ed25519
docker compose up -d --build
```

For a server build, the launcher selects the repository deploy key delivered
for that server and then runs Compose:

```bash
./up-server.sh
```

The server launcher loads its deploy key into a temporary SSH agent. Both SSH
identities are exposed only to the dependency-install build step. They are not
copied into the image or provided to the running container. A local developer
must have an SSH identity with read access to `ncats/IFX_Registry`;
`git ls-remote git@github.com:ncats/IFX_Registry.git` is a quick check.

At runtime, the base Compose file bind-mounts each YAML credential used by the
QA Browser read-only under `/app/src/use_cases/secrets`. Those YAML files remain
available even though `.dockerignore` excludes the secrets directory from the
build context. Mounting only the required YAMLs also keeps the server deploy key
out of the running container.
