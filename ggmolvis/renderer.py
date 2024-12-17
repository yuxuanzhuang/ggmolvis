import bpy
import tempfile
from IPython.display import display, Video, Image


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

    def render(self):
        # Render the scene
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

    def render(self):
        bpy.ops.render.render(animation=True)

    def display_in_notebook(self):
        return display(Video(url=self.filepath))
