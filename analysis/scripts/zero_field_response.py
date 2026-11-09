"""Read explicitly zero-field ORCA chi*T tensors, without using the scalar table."""
import re
import numpy as np

def zero_field_tensors(text):
    headers=list(re.finditer(r'Magnetic field \(Gauss\):\s*\(([^)]+)\)',text))
    found={}
    for index,header in enumerate(headers):
        field=[float(x) for x in header[1].split(',')]
        if not all(abs(x)<1e-12 for x in field):continue
        region=text[header.end():headers[index+1].start() if index+1<len(headers) else len(text)]
        for match in re.finditer(r'TEMPERATURE/K:\s+([\d.]+)\s*\nTensor in molecular frame \(cm3\*K/mol\):[^\n]*\n((?:[^\n]*\n){3})',region):
            temperature=float(match[1])
            matrix=np.array([[float(x) for x in line.split()] for line in match[2].splitlines()])
            assert matrix.shape==(3,3)
            assert temperature not in found, 'Ambiguous repeated zero-field temperature'
            found[temperature]=matrix
    assert found, 'No explicitly zero-field chi*T tensors found'
    return found
