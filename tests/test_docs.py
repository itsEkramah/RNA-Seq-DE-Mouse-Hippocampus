import re
from pathlib import Path

from scripts.data import ROOT


def test_readme_local_links_and_images_exist():
    readme = (ROOT / 'README.md').read_text(encoding='utf8')
    targets = re.findall(r'!?\[[^\]]*\]\(([^)]+)\)', readme)
    local = [target.split('#', 1)[0] for target in targets
             if not target.startswith(('https://', 'http://', 'mailto:'))]
    missing = [target for target in local if not (ROOT / target).exists()]
    assert not missing, f'Broken README links: {missing}'
