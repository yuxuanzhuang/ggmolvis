import bpy
import tempfile
from IPython.display import display, Video, Image
from contextlib import contextmanager, redirect_stdout, redirect_stderr
import os
import sys
from math import floor


@contextmanager
def BlenderProgressBar():
    """
    Context manager that shows a progress bar during Blender frame changes.
    Automatically removes the handler when exiting the context.

    Usage:
        with BlenderProgressBar():
            bpy.ops.render.render(animation=True)
    """

    def progress_callback(scene):
        # Get frame range
        start = scene.frame_start
        end = scene.frame_end
        current = scene.frame_current

        # Calculate progress
        total_frames = end - start + 1
        current_progress = current - start + 1
        percentage = (current_progress / total_frames) * 100

        # Create progress bar
        bar_width = 40
        filled_width = floor(bar_width * current_progress / total_frames)
        bar = "=" * filled_width + "-" * (bar_width - filled_width)

        # Create status string
        status = f"\rFrame {current}/{end} [{bar}] {percentage:.1f}%"

        # Print without newline and flush immediately
        print(status, end="", flush=True)

        # Print newline on last frame
        if current == end:
            print()

    try:
        # Register the handler
        bpy.app.handlers.frame_change_pre.append(progress_callback)
        yield
    finally:
        # Always remove the handler on exit
        if progress_callback in bpy.app.handlers.frame_change_pre:
            bpy.app.handlers.frame_change_pre.remove(progress_callback)


@contextmanager
def suppress_blender_output(match_string=None):
    """
    Context manager to suppress all output from Blender,
    including Python print statements. If a match_string
    is provided, lines containing the match_string will be printed.

    Usage:
        with suppress_blender_output(match_string="specific text"):
            bpy.ops.render.render()
    """
    # Save original file descriptors
    original_stdout_fd = os.dup(1)
    original_stderr_fd = os.dup(2)
    original_stdout = sys.stdout
    original_stderr = sys.stderr

    # Create a temporary file to capture the output
    temp_file = tempfile.TemporaryFile(mode="w+t")

    try:
        # Replace file descriptors at OS level
        os.dup2(temp_file.fileno(), 1)  # stdout
        os.dup2(temp_file.fileno(), 2)  # stderr

        # Redirect Python's stdout/stderr
        sys.stdout = temp_file
        sys.stderr = temp_file

        yield

    finally:
        # Restore original file descriptors
        os.dup2(original_stdout_fd, 1)
        os.dup2(original_stderr_fd, 2)

        # Restore Python's stdout/stderr
        sys.stdout = original_stdout
        sys.stderr = original_stderr

        # Clean up file descriptors
        os.close(original_stdout_fd)
        os.close(original_stderr_fd)

        # Read and print the captured output
        temp_file.seek(0)
        for line in temp_file:
            if match_string and match_string in line:
                print(line, end="")

        # Close the temporary file
        temp_file.close()


class Renderer:
    def __init__(
        self,
        resolution=(1280, 720),
        filepath=None,
        engine="EEVEE",
        cycles_device="GPU",
        format="PNG",
        samples: int = 64,
    ):
        self.resolution = resolution
        self.filepath = (
            filepath or tempfile.NamedTemporaryFile(suffix=f".{format}").name
        )
        self.engine = engine
        self.cycles_device = cycles_device
        self.samples = samples
        print("Rendering to:", self.filepath)
        self.format = format

    @property
    def cycles_device(self) -> str:
        return bpy.context.scene.cycles.device

    @cycles_device.setter
    def cycles_device(self, value: str) -> None:
        bpy.context.scene.cycles.device = value

    @property
    def engine(self) -> str:
        return bpy.context.scene.render.engine

    @engine.setter
    def engine(self, value: str) -> None:
        if value == "EEVEE":
            value = "BLENDER_EEVEE_NEXT"
        elif value == "WORKBENCH":
            value = "BLENDER_WORKBENCH"
        bpy.context.scene.render.engine = value

    @property
    def format(self):
        return bpy.context.scene.render.image_settings.file_format

    @format.setter
    def format(self, value):
        bpy.context.scene.render.image_settings.file_format = value

    @property
    def resolution(self):
        return (
            bpy.context.scene.render.resolution_x,
            bpy.context.scene.render.resolution_y,
        )

    @resolution.setter
    def resolution(self, value):
        bpy.context.scene.render.resolution_x = value[0]
        bpy.context.scene.render.resolution_y = value[1]

    @property
    def filepath(self):
        return bpy.context.scene.render.filepath

    @filepath.setter
    def filepath(self, value):
        bpy.context.scene.render.filepath = value

    @property
    def samples(self) -> int:
        return bpy.context.scene.eevee.taa_render_samples

    @samples.setter
    def samples(self, value: int) -> None:
        bpy.context.scene.eevee.taa_render_samples = value

    def render(self, verbose=False):
        if verbose:
            bpy.ops.render.render(write_still=True)
        else:
            with suppress_blender_output():
                bpy.ops.render.render(write_still=True)

    def display_in_notebook(self):
        # Load the image from the filepath
        return display(Image(self.filepath))


class MovieRenderer(Renderer):
    def __init__(
        self,
        resolution=(640, 360),
        filepath=None,
        engine="WORKBENCH",
        samples=64,
        cycles_device="GPU",
    ):
        super().__init__(
            resolution=resolution,
            filepath=filepath,
            format="FFMPEG",
            samples=samples,
            cycles_device=cycles_device,
            engine=engine,
        )
        print("Rendering to:", bpy.context.scene.render.filepath)
        self.video_format = "MPEG4"
        self.video_codec = "H264"
        self.video_quality = "HIGH"

    @property
    def video_forrmat(self) -> str:
        return bpy.context.scene.render.ffmpeg.format

    @video_forrmat.setter
    def video_format(self, value: str) -> None:
        bpy.context.scene.render.ffmpeg.format = value

    @property
    def video_codec(self) -> str:
        return bpy.context.scene.render.ffmpeg.codec

    @video_codec.setter
    def video_codec(self, value: str) -> None:
        bpy.context.scene.render.ffmpeg.codec = value

    @property
    def video_quality(self) -> str:
        return bpy.context.scene.render.ffmpeg.constant_rate_factor

    @video_quality.setter
    def video_quality(self, value: str) -> None:
        bpy.context.scene.render.ffmpeg.constant_rate_factor = value

    def render(self, verbose=False):
        if verbose:
            bpy.ops.render.render(animation=True)
        else:
            # First establish the progress bar context
            with BlenderProgressBar():
                # Then suppress only Blender's output inside
                with suppress_blender_output():
                    bpy.ops.render.render(animation=True)

    def display_in_notebook(self):
        return display(Video(url=self.filepath))
