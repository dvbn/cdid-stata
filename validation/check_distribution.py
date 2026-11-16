"""Check the Stata distribution using only the Python standard library."""
import hashlib
import json
import re
from pathlib import Path
from urllib.parse import unquote, urlsplit


ROOT = Path(__file__).resolve().parents[1]


def require(condition, message):
    if not condition:
        raise SystemExit(message)


def check_digest(path, expected):
    require(path.is_file() and not path.is_symlink(), f"Missing regular file: {path}")
    data = path.read_bytes()
    require(len(data) == expected["bytes"], f"Wrong byte count: {path}")
    require(hashlib.sha256(data).hexdigest() == expected["sha256"], f"Wrong SHA-256: {path}")


def check_catalog(directory, name):
    pkg = directory / f"{name}.pkg"
    toc = directory / "stata.toc"
    require(any(line.startswith(f"p {name} ") for line in toc.read_text().splitlines()),
            f"{name} is missing from {toc}")
    paths = set()
    for line in pkg.read_text().splitlines():
        fields = line.split()
        if not fields or fields[0] not in ("f", "F", "g", "G"):
            continue
        member = fields[-1]
        path = directory / member
        require(path.resolve().is_relative_to(directory) and path.is_file(),
                f"Invalid install catalog entry: {pkg}: {member}")
        paths.add(member)
    return paths


def check_links(paths):
    count = 0
    for path in paths:
        if path.suffix != ".md":
            continue
        content = re.sub(r"^```[^\n]*\n.*?^```\s*$", "", path.read_text(), flags=re.M | re.S)
        for match in re.finditer(r'\]\(([^\s)]+)\)|<img\b[^>]*\bsrc="([^"]+)"', content):
            link = urlsplit(next(value for value in match.groups() if value))
            if link.scheme or link.netloc or not link.path:
                continue
            target = (path.parent / unquote(link.path)).resolve()
            require(target.is_relative_to(ROOT) and target.exists(),
                    f"Broken local link: {path.relative_to(ROOT)}: {link.path}")
            count += 1
    return count


def main():
    manifest = json.loads((ROOT / "RELEASE_MANIFEST.json").read_text())
    expected = manifest["output_files"]
    paths = {p.relative_to(ROOT).as_posix() for p in ROOT.rglob("*")
             if (p.is_file() or p.is_symlink()) and ".git" not in p.parts}
    require(paths == set(expected) | {"RELEASE_MANIFEST.json"},
            f"Manifest file set differs: {sorted(paths ^ (set(expected) | {'RELEASE_MANIFEST.json'}))}")
    for name, digest in expected.items():
        check_digest(ROOT / name, digest)
    native = manifest["native_library"]
    check_digest(ROOT / native["file"], native)
    for name, digest in manifest["ordinary_dependency"]["files"].items():
        check_digest(ROOT / "csdid" / name, digest)
    install = check_catalog(ROOT, "cdid")
    dep = check_catalog(ROOT / "csdid", "csdid")
    require({"cdid.ado", "cdid.sthlp", native["file"]}.issubset(install),
            "The main installation is incomplete")
    require("pkg/csdid.ado" in dep, "The ordinary-mode installation is incomplete")
    version = manifest["package"]["version"]
    first = (ROOT / "cdid.ado").read_text().splitlines()[0].split()
    require(first[:3] == ["*!", "version", version], "Version differs from cdid.ado")
    links = check_links([ROOT / name for name in sorted(expected)])
    print(f"Verified {len(expected)} files, {len(install)} cdid catalog entries, "
          f"{len(dep)} csdid catalog entries and {links} local Markdown links.")


if __name__ == "__main__":
    main()
