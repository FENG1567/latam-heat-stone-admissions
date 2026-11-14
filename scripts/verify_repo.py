"""Verify package integrity and input availability without fitting a model."""
from pathlib import Path
import argparse,ast,csv,hashlib,json,re,sys

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--skip-manifest',action='store_true');args=parser.parse_args()
    root=Path(__file__).resolve().parents[1]
    needed=['README.md','README_中文.md','LICENSE','LICENSE-DATA.md','THIRD_PARTY_NOTICES.md','CITATION.cff','scripts/prepare_inputs.py','scripts/analysis.py','config/module_manifest.json','data/analysis_inputs/manifest.json','code/publication/run_publication.py','environment/requirements-publication.txt','paper/manuscript.md','paper/supplementary_tables.md']
    for name in needed:
        if not (root/name).is_file():raise FileNotFoundError(name)
    expected='284886236fa12d4ab8f6968bcb7f85244cc2f0337dbf07688e4b9e1181355b47'
    assert sha(root/'data/frozen/frozen_numbers.md')==expected,'Final numerical freeze changed'
    assert (root/'data/frozen/frozen_numbers.sha256').read_bytes()==(expected+'  frozen_numbers.md\n').encode()
    denied={'.cdsapirc','.env','id_rsa','id_ed25519','credentials.json'}
    patterns=[re.compile(r'gh[pousr]_[A-Za-z0-9]{30,}'),re.compile(r'AKIA[A-Z0-9]{16}'),re.compile(r'-----BEGIN (?:RSA |OPENSSH )?PRIVATE KEY-----')]
    py_count=0;files=[]
    for p in root.rglob('*'):
        if not p.is_file() or any(x in p.relative_to(root).parts for x in ['build','.venv','__pycache__','.git']):continue
        rel=p.relative_to(root).as_posix();files.append(p)
        if p.name in denied:raise RuntimeError('Credential-like file must not be published: '+rel)
        if p.stat().st_size>25*1024*1024:raise RuntimeError('File exceeds browser upload limit: '+rel)
        if p.suffix=='.py':ast.parse(p.read_text(encoding='utf-8-sig'),filename=rel);py_count+=1
        if p.suffix in ['.py','.R','.json','.sh','.md','.txt','.yml','.cff']:
            text=p.read_text(encoding='utf-8-sig')
            for pat in patterns:
                if pat.search(text):raise RuntimeError('Credential signature found in '+rel)
    spec=json.loads((root/'data/analysis_inputs/manifest.json').read_text(encoding='utf-8'))
    for item in spec['inputs']:
        h=hashlib.sha256();size=0
        for part in item['parts']:
            p=(root/part['path']).resolve();assert p.is_relative_to(root)
            assert p.is_file() and p.stat().st_size==part['bytes'] and sha(p)==part['sha256'],part['path']
            with p.open('rb') as f:
                for b in iter(lambda:f.read(1024*1024),b''):h.update(b);size+=len(b)
        assert size==item['bytes'] and h.hexdigest()==item['sha256'],item['target']
    csv_files=list((root/'paper/source_data').glob('*.csv'));assert len(csv_files)==87
    index=list(csv.DictReader((root/'data/source_data_index.csv').open(encoding='utf-8')));assert len(index)==87
    for row in index:assert sha(root/row['file'])==row['sha256']
    for name in ['manuscript.md','manuscript_中文.md','tables.md','tables_中文.md','supplementary_methods.md','supplementary_tables.md']:
        p=root/'paper'/name;s=p.read_text(encoding='utf-8')
        for link in re.findall(r'\]\(([^)]+)\)',s):
            if not re.match(r'https?://|mailto:',link):assert (p.parent/link).is_file(),(name,link)
    manifest=root/'MANIFEST.sha256'
    if manifest.exists() and not args.skip_manifest:
        for line in manifest.read_text(encoding='utf-8').splitlines():
            expected,rel=line.split('  ',1);p=root/rel
            assert p.is_file() and sha(p)==expected,'Package manifest mismatch: '+rel
    print(json.dumps({'status':'passed','source_data_csvs':87,'aggregate_inputs':len(spec['inputs']),'python_sources_parsed':py_count,'files_checked':len(files),'max_file_bytes':max(p.stat().st_size for p in files),'credential_signatures':0,'freeze_verified':True,'scope':'No statistical model refits'},indent=2))
if __name__=='__main__':main()
