version 1.0

workflow ab1_alignment_pipe {
    input{
        File reference_fasta
        String path_to_ab1_files
    }
    

    call run_tracy_alignment {
        input:
            ab1_files = ab1_files,
            reference_fasta = reference_fasta
    }
    
    call plot_alignment_general {
        input:
            out_aln_json = run_tracy_alignment.out_aln_json
            ref_basename = run_tracy_alignment.ref_basename
    }
    
    call plot_alignment_sabre {
        # input:
            
    }
}


task run_tracy_alignment {
    input{
        String ab1_files
        File reference_fasta 
    }

    command <<<
        tracy assemble -r ~{reference_fasta} ~{ab1_files} -o out.align
    >>>

    output {
        File out_aln_fa = 'out.align.fa'
        File out_aln_json = 'out.align.json'
    }
}


task plot_alignment_general {
    input{
        # File out_aln_fa
        File out_aln_json  
        String pdf_name = "out.align.general_plot.pdf"
        }
    
    command <<<
        python <<CODE
            import json
            import matplotlib.pyplot as plt
            import numpy as np
            import matplotlib.colors as mcolors


            with open("~{out_aln_json}") as f:
                data = json.load(f)

            msa = data["msa"]

            seqs = [s["align"] for s in msa]
            labels = [s["traceFileName"] if s["traceFileName"] else "REFERENCE" for s in msa]

            clear_aln = []
            for s in msa: 
                clear_aln.append(' '*int(s['leadingGaps']) + s['align'] + ' '*int(s['trailingGaps']))

            alphabet = [' ', '-', 'A', 'T', 'G', 'C', 'N']
            char_to_int = {char: i for i, char in enumerate(alphabet)}

            array_aln = np.array([[char_to_int.get(base, 0) for base in seq] for seq in clear_aln])

            xlim_left = min([int(seq['leadingGaps']) for seq in msa[:-1]]) - 50
            xlim_right = max([len(msa[-1]['align']) - int(seq['trailingGaps']) for seq in msa[:-1]]) + 50


            colors = [
                "white",            # 0 -- ' ' (padding)
                "pink",             # 1 -- without coverage
                "lightgreen",       # 2 -- with coverage
                "black"             # 3 -- mismatch
            ]

            color_labels = [
                'padding',
                'no coverage',
                'covered',
                'mismatch'
            ]

            cmap = mcolors.ListedColormap(colors)
            bounds = np.arange(len(colors) + 1) - 0.5 
            norm = mcolors.BoundaryNorm(bounds, cmap.N)

            plt.figure(figsize=(18, len(seqs) * 0.3))
            plt.xlim(xlim_left, xlim_right)

            # fill white
            ref_seq = array_aln[-1]
            visual_array = np.full_like(array_aln, fill_value=0)


            for j in range(array_aln.shape[1]):
                ref_base = array_aln[-1, j]

                has_coverage = any(
                    array_aln[i, j] != char_to_int[' ']
                    for i in range(len(array_aln) - 1)  # exclude reference
                )

                has_mismatch = any(
                    (array_aln[i, j] not in [char_to_int[' ']]) and
                    (array_aln[i, j] != ref_base)
                    for i in range(len(array_aln) - 1)
                )

                for i in range(array_aln.shape[0]):
                    current_base = array_aln[i, j]
                    # looking at the reference
                    if i == len(array_aln) - 1:
                        # reference + no coverage -- PINK
                        if has_coverage == False:
                            visual_array[i, j] = 1
                        # reference + coverage + match -- GREEN
                        elif has_mismatch == False:
                            visual_array[i, j] = 2
                        # reference + coverage + mismatches -- BLACK
                        else:
                            visual_array[i, j] = 3
                    
                    # looking at the sequences:
                    else:
                        # seq + no coverage -- WHITE
                        if not has_coverage:
                            visual_array[i, j] = 0
                        else:
                            # seq + coverage + padding -- WHITE
                            if current_base == char_to_int[' ']:
                                continue                
                            else:
                                # seq + coverage + match -- GREEN
                                if current_base == ref_base:
                                    visual_array[i, j] = 2
                                # seq + coverage + mismatch -- BLACK
                                else:
                                    visual_array[i, j] = 3
                            
            plt.imshow(visual_array, cmap=cmap, norm=norm, aspect='auto', interpolation='none')

            plt.yticks(ticks=range(len(labels)), labels=labels)
            plt.xlabel("Position in alignment")
            plt.colorbar(
                ticks=list(range(len(color_labels))),
                format=plt.FuncFormatter(lambda val, _: color_labels[int(val)]),
            )
            plt.title("Alignment visualization")
            plt.tight_layout()
            plt.savefig("~{pdf_name}", dpi=500)
            # plt.show()
        CODE
        
    >>>
    
    output{
        File pdf_plot = pdf_name 
    }
}

task plot_alignment_sabre {
    
}