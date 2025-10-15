import os
import subprocess

def generate_tikz_file(result, filename, parameter_name, parameter_value, tikz_dir):
    """Generate TikZ visualization from results."""
    
    if not result or not result['success']:
        print(f"Cannot generate TikZ file for {filename}: optimization failed")
        return
    
    epsilons = result['epsilons']
    weights = result['weights']
    
    edge_mapping = [
        (1, 2), (1, 4), (1, 7), (1, 10),
        (2, 3), (2, 4),
        (3, 4), (3, 5),
        (4, 5),
        (5, 6),
        (6, 7), (6, 9),
        (7, 8),
        (8, 9), (8, 10),
        (9, 10)
    ]
    
    weight_names = [
        "one", "two", "three", "four", "five", "six", "seven", "eight", 
        "nine", "ten", "eleven", "twelve", "thirteen", "fourteen", "fifteen", "sixteen"
    ]
    
    privacy_names = [
        "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten"
    ]
    
    pastel_colors = [
        "pastel1", "pastel2", "pastel3", "pastel4", "pastel5",
        "pastel6", "pastel7", "pastel8", "pastel9", "pastel10"
    ]
    
    tikz_content = f"""\\documentclass[border=5mm]{{standalone}}
\\usepackage{{tikz}}
\\usepackage{{xcolor}}
\\usetikzlibrary{{positioning}}

\\definecolor{{pastel1}}{{RGB}}{{255,182,193}}
\\definecolor{{pastel2}}{{RGB}}{{173,216,230}}
\\definecolor{{pastel3}}{{RGB}}{{144,238,144}}
\\definecolor{{pastel4}}{{RGB}}{{255,218,185}}
\\definecolor{{pastel5}}{{RGB}}{{221,160,221}}
\\definecolor{{pastel6}}{{RGB}}{{255,255,224}}
\\definecolor{{pastel7}}{{RGB}}{{176,224,230}}
\\definecolor{{pastel8}}{{RGB}}{{250,128,114}}
\\definecolor{{pastel9}}{{RGB}}{{152,251,152}}
\\definecolor{{pastel10}}{{RGB}}{{255,160,122}}

\\begin{{document}}
\\begin{{tikzpicture}}[
    node distance=5cm
]

% Privacy parameters for {parameter_name} = {parameter_value}
"""

    for i, privacy_name in enumerate(privacy_names):
        tikz_content += f"\\newcommand{{\\privacy{privacy_name}}}{{{epsilons[i]:.3f}}}\n"
    
    tikz_content += f"""
% Edge weights for {parameter_name} = {parameter_value}
"""

    for i, weight_name in enumerate(weight_names):
        if i < len(weights):
            edge_info = edge_mapping[i] if i < len(edge_mapping) else (i+1, i+2)
            tikz_content += f"\\newcommand{{\\weight{weight_name}}}{{{weights[i]:.3f}}}\n"
    
    tikz_content += """
% Node sizes based on privacy params
"""

    for i, privacy_name in enumerate(privacy_names):
        tikz_content += f"\\pgfmathsetmacro{{\\size{privacy_name}}}{{0.6 + 5 * \\privacy{privacy_name}}}\n"
    
    tikz_content += """
% Edge thicknesses
"""

    for i, weight_name in enumerate(weight_names):
        if i < len(weights):
            if i == 5:
                tikz_content += f"\\pgfmathsetmacro{{\\thick{weight_name}}}{{max(0.1, 0.3 + 0.5 * \\weight{weight_name} * 10)}}\n"
            else:
                tikz_content += f"\\pgfmathsetmacro{{\\thick{weight_name}}}{{0.3 + 0.5 * \\weight{weight_name} * 10}}\n"
    
    tikz_content += """
% Nodes
"""

    angles = [0, 36, 72, 108, 144, 180, 216, 252, 288, 324]
    for i, (privacy_name, pastel_color, angle) in enumerate(zip(privacy_names, pastel_colors, angles)):
        node_num = i + 1
        tikz_content += f"\\node[circle, draw=black, fill={pastel_color}, minimum size=\\size{privacy_name} cm, font=\\Large\\bfseries, line width=1pt] ({node_num}) at ({angle}:5) {{{node_num}}};\n"
    
    tikz_content += """
% Edges
"""

    for i, (weight_name, edge) in enumerate(zip(weight_names, edge_mapping)):
        if i >= len(weights):
            break
            
        node1, node2 = edge
        
        if weights[i] < 0.004:
            continue
        else:
            tikz_content += f"\\draw[line width=\\thick{weight_name} pt, draw=gray] ({node1}) -- ({node2});\n"
            
            if edge in [(3, 5), (6, 9), (8, 10)]:
                if edge == (3, 5):
                    position = "pos=0.55, below left"
                elif edge == (6, 9):
                    position = "pos=0.75"
                else:
                    position = "pos=0.6, above"
            else:
                position = "midway"
            
            tikz_content += f"\\path ({node1}) -- ({node2}) node[{position}, fill=white, inner sep=2pt, font=\\normalsize\\bfseries] {{\\pgfmathprintnumber[fixed, precision=3]{{\\weight{weight_name}}}}};\n"
    
    tikz_content += """
\\end{tikzpicture}
\\end{document}
"""

    os.makedirs(tikz_dir, exist_ok=True)
    
    with open(filename, 'w', encoding='utf-8') as f:
        f.write(tikz_content)
    
    print(f"Generated TikZ file: {filename}")


def generate_all_tikz_files(results, parameter_name, file_prefix, tikz_dir):
    """Generate TikZ files for all successful runs."""
    
    print(f"\nGenerating TikZ files for {parameter_name} comparison...")
    
    for param_value, result in results.items():
        if result and result['success']:
            if isinstance(param_value, float):
                param_str = f"{param_value:.1f}".replace('.', '_')
            else:
                param_str = str(param_value)
            
            filename = f"{tikz_dir}/{file_prefix}_{param_str}.tex"
            generate_tikz_file(result, filename, parameter_name, param_value, tikz_dir)


def generate_tikz_pdf(tex_filename):
    """Try to compile TikZ to PDF with pdflatex."""
    try:
        tex_dir = os.path.dirname(tex_filename)
        tex_base = os.path.basename(tex_filename)
        
        result = subprocess.run(['pdflatex', '-interaction=nonstopmode', tex_base], 
                              cwd=tex_dir, capture_output=True, text=True, timeout=30)
        
        if result.returncode == 0:
            pdf_filename = tex_filename.replace('.tex', '.pdf')
            print(f"  -> PDF generated: {pdf_filename}")
            
            base_name = tex_base.replace('.tex', '')
            for ext in ['.aux', '.log', '.synctex.gz', '.fls', '.fdb_latexmk']:
                aux_file = os.path.join(tex_dir, base_name + ext)
                if os.path.exists(aux_file):
                    try:
                        os.remove(aux_file)
                    except:
                        pass
        else:
            print(f"  -> PDF compilation failed for {tex_filename}")
            
    except FileNotFoundError:
        print(f"  -> pdflatex not found, skipping PDF generation")
    except Exception as e:
        print(f"  -> Error generating PDF: {e}")
