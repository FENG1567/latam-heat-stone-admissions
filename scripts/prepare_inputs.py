"""Restore checksum-bound aggregate inputs; never refit or change data values."""
from pathlib import Path
import argparse,hashlib,json,os,tempfile

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def within(root,relative):
    path=(root/relative).resolve()
    if not path.is_relative_to(root.resolve()):raise ValueError('Path escapes target root: '+relative)
    return path

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest',type=Path)
    parser.add_argument('--target-root','--outdir',dest='target',type=Path)
    parser.add_argument('--verify-only',action='store_true')
    args=parser.parse_args();repo=Path(__file__).resolve().parents[1]
    manifest=args.manifest or repo/'data/analysis_inputs/manifest.json'
    target=args.target or repo/'build/analysis_workspace'
    spec=json.loads(manifest.read_text(encoding='utf-8'));records=[]
    for item in spec['inputs']:
        for part in item['parts']:
            path=within(repo,part['path'])
            if not path.is_file() or path.stat().st_size!=part['bytes'] or sha(path)!=part['sha256']:
                raise ValueError('Input part missing or changed: '+part['path'])
        destination=within(target,item['target'])
        if args.verify_only:records.append({'target':item['target'],'status':'parts_verified'});continue
        destination.parent.mkdir(parents=True,exist_ok=True)
        if destination.exists():
            if sha(destination)!=item['sha256']:raise ValueError('Existing input differs; keep it safe and use a fresh target root: '+str(destination))
            records.append({'target':item['target'],'status':'already_verified'});continue
        fd,tmp=tempfile.mkstemp(prefix='.input_',dir=destination.parent);temp=Path(tmp)
        try:
            with os.fdopen(fd,'wb') as out:
                for part in item['parts']:
                    with within(repo,part['path']).open('rb') as src:
                        for block in iter(lambda:src.read(1024*1024),b''):out.write(block)
            if temp.stat().st_size!=item['bytes'] or sha(temp)!=item['sha256']:raise ValueError('Restored input checksum mismatch: '+item['target'])
            os.replace(temp,destination)
        finally:temp.unlink(missing_ok=True)
        records.append({'target':item['target'],'status':'restored_byte_identically'})
    print(json.dumps({'status':'passed','inputs':len(records),'mode':'verify_parts' if args.verify_only else 'restore_original_bytes','target_root':str(target)},indent=2))
if __name__=='__main__':main()
