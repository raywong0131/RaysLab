"""Save a supplied scientific figure without changing its numerical display."""
from .output_paths import artifact_path


def save_figure_formats(figure, root, stem, *, formats=("png", "pdf"), pdf_group=None,
                        dpi=None, **save_options):
    paths = []
    for extension in formats:
        path = artifact_path(root, f"{stem}.{extension}",
                             group=pdf_group if extension == "pdf" else None)
        options = dict(save_options)
        if dpi is not None:
            options["dpi"] = dpi
        figure.savefig(path, **options)
        paths.append(path)
    return paths
